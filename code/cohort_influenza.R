#' Get flu diagnosed patients, from scrach, from condition occurence table.
#' Does not take an existing cohort as a parameter
#'
#' @param min_date 
#' @param max_date 
#'
#' @return Dataframe with person ids and earliest flu diagnosis date
#' @export
#'
#' @examples
#' my_flu_cohort <- get_flu_evidence(as.Date("2019-01-01"), as.Date("2020-01-01))
#' 
get_flu_evidence <- function(cohort_tbl, min_date, max_date) {
  
  flu_conditions <-
    load_codeset("dx_influenza") %>% 
    left_join(cdm_tbl("condition_occurrence") %>% 
                select(condition_source_concept_id, person_id, site, condition_start_date) %>% 
                filter(condition_start_date < max_date, condition_start_date >= min_date) %>% 
                inner_join(cohort_tbl %>% select(person_id), by="person_id"), 
              by=c("concept_id"="condition_source_concept_id")) %>% 
    compute_new(indexes=c("person_id"))
  
  flu_tests <- 
    load_codeset("lab_influenza_filtered") %>% 
    left_join(cdm_tbl("measurement_labs") %>% 
                select(measurement_concept_id, person_id, site, measurement_date, value_as_concept_id) %>% 
                filter(measurement_date < max_date, measurement_date >= min_date) %>% 
                inner_join(cohort_tbl %>% select(person_id), by="person_id"), 
              by=c("concept_id"="measurement_concept_id")) %>% 
    compute_new(indexes=c("person_id"))
  
  # Should be a full join
  flu_tests %>% 
    filter(value_as_concept_id == 9191) %>% 
    group_by(person_id) %>% # This part just gets the patient's earliest flu evidence, but instead, could just mark a new column as "first" that could later be filtered
    slice_min(measurement_date, with_ties = FALSE) %>% 
    ungroup() %>% 
    select(person_id, earliest_flu_test = measurement_date) %>% 
    compute_new(indexes=c("person_id")) %>% 
    full_join(flu_conditions %>% 
                group_by(person_id) %>% 
                slice_min(condition_start_date, with_ties = FALSE) %>% # Same here, instead of slicing the min, could just identify the min
                ungroup() %>% 
                select(person_id, flu_dx_date = condition_start_date) %>% 
                compute_new(indexes=c("person_id")), by=c("person_id")) %>% 
    mutate(earliest_flu_evidence = case_when(is.na(earliest_flu_test) ~ flu_dx_date,
                                             is.na(flu_dx_date) ~ earliest_flu_test,
                                             flu_dx_date <= earliest_flu_test ~ flu_dx_date,
                                             earliest_flu_test < flu_dx_date ~ earliest_flu_test)) %>% 
    return()
}

get_respiratory_evidence <- function(cohort_tbl, min_date, max_date) {
  # TODO Use new codeset
  noncovid_resp_conditions <-
    load_codeset("dx_noncov_resp_infections") %>% 
    inner_join(cohort %>% 
                 select(person_id) %>% 
                 left_join(cdm_tbl("condition_occurrence") %>% 
                             select(visit_occurrence_id, condition_occurrence_id, condition_source_concept_id, person_id, site, condition_start_date) %>% 
                             filter(condition_start_date < max_date, condition_start_date >= min_date),
                           by="person_id") %>% compute_new(indexes=c("person_id", "condition_source_concept_id")), 
               by=c("concept_id"="condition_source_concept_id")) %>% 
    compute_new(indexes=c("person_id")) %>% 
    group_by(person_id) %>% 
    slice_min(condition_start_date) %>% 
    ungroup() %>% 
    compute_new(indexes=c("person_id"))
  
  cohort_tbl %>% 
    left_join(noncovid_resp_conditions %>% 
                select(person_id, earliest_resp_evidence = condition_start_date),
              by="person_id") %>% 
    return()
  
  
}


build_ili_controls <- function(cohort_antijoin) {
  ili_codes <-
    load_codeset('ili_conditions')
  
  conditions_dates <-
    cdm_tbl("condition_occurrence") %>% 
    filter(condition_start_date > as.Date("2020-03-01"), condition_start_date < as.Date("2023-09-01")) %>% 
    select(person_id, condition_start_date, condition_concept_id, condition_occurrence_id)
  
  ili_conditions <- ili_codes %>% 
    left_join(conditions_dates, by=c("concept_id"="condition_concept_id")) %>% 
    compute_new(indexes=c("person_id"))
  
  patients_first_date <-
    ili_conditions %>% 
    group_by(person_id) %>% 
    slice_min(condition_start_date, with_ties = FALSE) %>% 
    ungroup() %>% 
    compute_new(indexes=c("person_id"))
  
  ili_pats_filtered <-
    patients_first_date %>% 
    mutate(index_date = condition_start_date) %>% 
    flag_visits_prior(cdm_tbl("visit_occurrence")) %>%
    # flag_visit_later(cdm_tbl("visit_occurrence")) %>% 
    filter(visit_within_2_years == TRUE) %>% 
    compute_new(indexes=c("person_id"))
  
  
  ili_pats_antijoined <-
    ili_pats_filtered %>% 
    anti_join(cohort_antijoin %>% select(person_id),
              by="person_id") %>% 
    compute_new(indexes=c("person_id"))
  
  ili_pats_antijoined %>% 
    return()
}
