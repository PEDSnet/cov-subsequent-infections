#' Get flu diagnosed patients for a cohort, returns earliest flu evidence between dates given
#' Cohort can just be some subset of the person_tbl, or the whole thing if necessary for finding patients
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
get_flu_evidence_NN <- function(cohort_tbl, 
                                min_date, max_date, 
                                test_only = FALSE,
                                flu_lab_codeset = "lab_influenza_complete", 
                                cohort_1_label) {
  
  flu_conditions <-
    load_codeset("dx_influenza_final") %>% 
    left_join(cdm_tbl("condition_occurrence") %>% 
                select(condition_source_concept_id, person_id, site, condition_start_date) %>% 
                filter(condition_start_date < max_date, condition_start_date >= min_date) %>% 
                inner_join(cohort_tbl %>% select(person_id), by="person_id"), 
              by=c("concept_id"="condition_source_concept_id")) %>% 
    compute_new(indexes=c("person_id"))
  flu_conditions %>% output_tbl(paste0(cohort_1_label, "_flu_dx"))
  
  flu_tests <- 
    load_codeset(flu_lab_codeset) %>% 
    left_join(cdm_tbl("measurement_labs") %>% 
                select(measurement_concept_id, person_id, site, measurement_date, value_as_concept_id, value_source_value) %>% 
                filter(measurement_date < max_date, measurement_date >= min_date) %>% 
                inner_join(cohort_tbl %>% select(person_id), by="person_id"), 
              by=c("concept_id"="measurement_concept_id")) %>% 
    compute_new(indexes=c("person_id"))
  
  flu_tests %>% output_tbl(paste0(cohort_1_label, "_flu_lab_complete")) 
  
  if (test_only == TRUE) {
    flu_tests %>% 
      filter(!grepl("undetected", value_source_value, ignore.case = TRUE)) %>%
      filter(!grepl("non-detected", value_source_value, ignore.case = TRUE)) %>%
      filter(!grepl("not detected", value_source_value, ignore.case = TRUE)) %>%
      filter(grepl("detected", value_source_value, ignore.case = TRUE)|
               grepl("positive", value_source_value, ignore.case = TRUE)|
               grepl("RNA detected", value_source_value, ignore.case = TRUE)|
               grepl("pos", value_source_value, ignore.case = TRUE)|
               value_as_concept_id %in% c(9191L, 2000001526)) %>%
      # !grepl("possible", value_source_value, ignore.case = TRUE)) %>%
      # filter(!grepl("faint pos", value_source_value, ignore.case = TRUE)) %>%
      # filter(!grepl("faintly pos", value_source_value, ignore.case = TRUE)) %>%
      filter(!grepl("negative", value_source_value, ignore.case = TRUE)) %>%
      filter(!grepl("presumptive", value_source_value, ignore.case = TRUE)) %>%
      # filter(!grepl("weakly pos", value_source_value, ignore.case = TRUE)) %>%
      # filter(!grepl("weakly positive", value_source_value, ignore.case = TRUE)) %>%
      filter(!grepl("SARS-CoV-2-Detected", value_source_value, ignore.case = TRUE)) %>%
      group_by(person_id) %>%
      slice_min(measurement_date, with_ties = FALSE) %>% 
      ungroup() %>% 
      select(person_id, earliest_flu_test = measurement_date,
             flu_lab_concept_id = concept_id,
             flu_value_concept_id = value_as_concept_id) %>% 
      compute_new(indexes=c("person_id")) %>% 
      mutate(earliest_flu_evidence = earliest_flu_test) %>% 
      return()
  } else{
    flu_tests %>% 
      filter(!grepl("undetected", value_source_value, ignore.case = TRUE)) %>%
      filter(!grepl("non-detected", value_source_value, ignore.case = TRUE)) %>%
      filter(!grepl("not detected", value_source_value, ignore.case = TRUE)) %>%
      filter(grepl("detected", value_source_value, ignore.case = TRUE)|
               grepl("positive", value_source_value, ignore.case = TRUE)|
               grepl("RNA detected", value_source_value, ignore.case = TRUE)|
               grepl("pos", value_source_value, ignore.case = TRUE)|
               value_as_concept_id %in% c(9191L, 2000001526)) %>%
      # !grepl("possible", value_source_value, ignore.case = TRUE)) %>%
      # filter(!grepl("faint pos", value_source_value, ignore.case = TRUE)) %>%
      # filter(!grepl("faintly pos", value_source_value, ignore.case = TRUE)) %>%
      filter(!grepl("negative", value_source_value, ignore.case = TRUE)) %>%
      filter(!grepl("presumptive", value_source_value, ignore.case = TRUE)) %>%
      # filter(!grepl("weakly pos", value_source_value, ignore.case = TRUE)) %>%
      # filter(!grepl("weakly positive", value_source_value, ignore.case = TRUE)) %>%
      filter(!grepl("SARS-CoV-2-Detected", value_source_value, ignore.case = TRUE)) %>%
      group_by(person_id) %>% 
      slice_min(measurement_date, with_ties = FALSE) %>% 
      ungroup() %>% 
      select(person_id, earliest_flu_test = measurement_date) %>% 
      compute_new(indexes=c("person_id")) %>% 
      full_join(flu_conditions %>% 
                  group_by(person_id) %>% 
                  slice_min(condition_start_date, with_ties = FALSE) %>% 
                  ungroup() %>% 
                  select(person_id, flu_dx_date = condition_start_date) %>% 
                  compute_new(indexes=c("person_id")), by=c("person_id")) %>% 
      mutate(earliest_flu_evidence = case_when(is.na(earliest_flu_test) ~ flu_dx_date, # only dx, no positive test
                                               is.na(flu_dx_date) ~ earliest_flu_test, # 
                                               flu_dx_date <= earliest_flu_test ~ flu_dx_date,
                                               earliest_flu_test < flu_dx_date ~ earliest_flu_test)) %>% 
      return()
  }
  
  
}

get_respiratory_evidence <- function(cohort_tbl, min_date, max_date) {
  noncovid_resp_conditions <-
    load_codeset("dx_other_resp_final") %>% 
    inner_join(cohort_tbl %>% 
                 select(person_id) %>% 
                 left_join(cdm_tbl("condition_occurrence") %>% 
                             select(visit_occurrence_id, condition_occurrence_id, condition_source_concept_id, person_id, site, condition_start_date) %>% 
                             filter(condition_start_date < max_date, condition_start_date >= min_date),
                           by="person_id") %>% compute_new(indexes=c("person_id", "condition_source_concept_id")), 
               by=c("concept_id"="condition_source_concept_id")) %>% 
    compute_new(indexes=c("person_id")) %>% 
    group_by(person_id) %>% 
    slice_min(condition_start_date, with_ties = FALSE) %>% 
    ungroup() %>% 
    compute_new(indexes=c("person_id"))
  
  noncovid_resp_conditions %>% 
    return()
  # 
  # cohort_tbl %>% 
  #   left_join(noncovid_resp_conditions %>% 
  #               select(person_id, earliest_resp_evidence = condition_start_date),
  #             by="person_id") %>% 
  #   return()
  
  
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

get_other_infection_labs <- function(cohort_tbl, min_date, max_date, cohort_1_label) {
  
  # flu_test <- cohort_tbl %>% 
  #   filter(ce_event ==) #filter out flu cohort, test without flu diagnosis
  
  # if that patient has a flu diagnosis within 14 days
  flu_dx <- results_tbl(paste0(cohort_1_label, "_flu_dx"))
  test_and_dx <- flu_test %>% 
    inner_join(flu_dx, by="person_id") %>%
    filter(ce_date <= condition_start_date + days(14), 
           ce_date >= condition_start_date - days(14)) %>%
    distinct(person_id) %>% 
    mutate(pos_and_dx = TRUE)
  
  test_no_dx <- flu_test %>% 
    anti_join(test_and_dx %>% select(person_id), by="person_id") %>%
    distinct(person_id, ce_date) %>%
    mutate(pos_no_dx = TRUE)
  
  # if this person has another positive tests within 14 days
  test_no_dx_other_test_pos <- test_no_dx %>% 
    inner_join(cdm_tbl("measurement_labs") %>% 
                 select(measurement_concept_id, person_id, site, measurement_date, value_as_concept_id) %>% 
                 filter(measurement_date <= ce_date + days(14), 
                        measurement_date >= ce_date - days(14)), by = "person_id") %>%
    inner_join(load_codeset("other_respiratory_infections_lab"), by=c("measurement_concept_id" = "concept_id")) %>% 
    filter(!grepl("undetected", value_source_value, ignore.case = TRUE)) %>%
    filter(!grepl("non-detected", value_source_value, ignore.case = TRUE)) %>%
    filter(!grepl("not detected", value_source_value, ignore.case = TRUE)) %>%
    filter(grepl("detected", value_source_value, ignore.case = TRUE)|
             grepl("positive", value_source_value, ignore.case = TRUE)|
             grepl("RNA detected", value_source_value, ignore.case = TRUE)|
             grepl("pos", value_source_value, ignore.case = TRUE)|
             value_as_concept_id %in% c(9191L, 2000001526)) %>%
    # !grepl("possible", value_source_value, ignore.case = TRUE)) %>%
    # filter(!grepl("faint pos", value_source_value, ignore.case = TRUE)) %>%
    # filter(!grepl("faintly pos", value_source_value, ignore.case = TRUE)) %>%
    filter(!grepl("negative", value_source_value, ignore.case = TRUE)) %>%
    filter(!grepl("presumptive", value_source_value, ignore.case = TRUE)) %>%
    # filter(!grepl("weakly pos", value_source_value, ignore.case = TRUE)) %>%
    # filter(!grepl("weakly positive", value_source_value, ignore.case = TRUE)) %>%
    filter(!grepl("SARS-CoV-2-Detected", value_source_value, ignore.case = TRUE)) %>%
    distinct(person_id) %>%
    mutate(pos_no_dx_other_infections = TRUE) %>%
    compute_new(indexes=c("person_id"))
  
  final_cohort <- cohort_tbl %>% 
    left_join(test_and_dx, by="person_id") %>%
    left_join(test_no_dx, by="person_id") %>%
    left_join(test_no_dx_other_test_pos, by="person_id") %>%
    mutate(flu_status = case_when(test_no_dx_other_test_pos ~ "test_no_dx_other_infections_pos",
                                  test_no_dx & !test_no_dx_other_test_pos ~ "test_only",
                                  test_and_dx & !test_no_dx ~ "test_and_dx",
                                  ce_event ==  ~ "dx_only",
                                  TRUE ~ NA)) %>% return()
}
