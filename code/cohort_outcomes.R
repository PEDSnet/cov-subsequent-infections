### All infection code

#' Title
#'
#' @param cohort 
#' @param outcome_start_date 
#' @param outcome_end_date 
#'
#' @return
#' @export
#'
#' @examples
flag_outcome_infections <- function(cohort, outcome_start_date, outcome_end_date) {
  
  any_infection_conditions <-
    load_codeset("dx_noncov_any_infection") %>% 
    inner_join(cohort %>%  # shouldn't this be an inner join??
                select(person_id) %>% 
                left_join(cdm_tbl("condition_occurrence") %>% 
                            select(visit_occurrence_id, condition_occurrence_id, condition_source_concept_id, person_id, site, condition_start_date) %>% 
                            filter(condition_start_date < outcome_end_date, condition_start_date >= outcome_start_date),
                          by="person_id") %>% compute_new(indexes=c("person_id", "condition_source_concept_id")), 
              by=c("concept_id"="condition_source_concept_id")) %>% 
    compute_new(indexes=c("person_id"))
  
  ## Returning a list of ALL infections, not just tidy / first one
  any_infection_conditions %>% 
    return()
}


### Respiratory infection code
#' Title
#'
#' @param cohort 
#' @param outcome_start_date 
#' @param outcome_end_date 
#'
#' @return
#' @export
#'
#' @examples
flag_outcome_resp_infections <- function(cohort, outcome_start_date, outcome_end_date) {
  # TODO make an option in here that also includes covid infections but can also just flag later
  noncovid_resp_conditions <-
    load_codeset("dx_noncov_resp_infections") %>% 
    inner_join(cohort %>% 
                select(person_id) %>% 
                left_join(cdm_tbl("condition_occurrence") %>% 
                            select(visit_occurrence_id, condition_occurrence_id, condition_source_concept_id, person_id, site, condition_start_date) %>% 
                            filter(condition_start_date < outcome_end_date, condition_start_date >= outcome_start_date),
                          by="person_id") %>% compute_new(indexes=c("person_id", "condition_source_concept_id")), 
              by=c("concept_id"="condition_source_concept_id")) %>% 
    compute_new(indexes=c("person_id"))
  
  ## Returning a list of ALL infections, not just tidy / first one
  noncovid_resp_conditions %>% 
    return()
}

#' Title
#'
#' @param cohort 
#' @param outcome_window_start 
#' @param outcome_window_end 
#'
#' @return
#' @export
#'
#' @examples
flag_rsv_outcome <- function(cohort, outcome_start_date, outcome_end_date, test_only = FALSE) {
  
  rsv_conditions <-
    load_codeset("dx_rsv_final") %>% 
    left_join(cohort %>% 
                select(person_id) %>% 
                left_join(cdm_tbl("condition_occurrence") %>% 
                            select(visit_occurrence_id, condition_occurrence_id, condition_source_concept_id, person_id, site, condition_start_date) %>% 
                            filter(condition_start_date < outcome_end_date, condition_start_date >= outcome_start_date),
                          by="person_id") %>% compute_new(indexes=c("person_id", "condition_source_concept_id")), 
              by=c("concept_id"="condition_source_concept_id")) %>% 
    compute_new(indexes=c("person_id"))
  
  rsv_tests <- 
    load_codeset("lab_rsv_final") %>% 
    left_join(cohort %>% select(person_id) %>% 
                left_join(cdm_tbl("measurement_labs") %>% 
                            select(measurement_concept_id, person_id, site, measurement_date, value_as_concept_id) %>% 
                            filter(measurement_date < outcome_end_date, measurement_date >= outcome_start_date),
                          by="person_id"), 
              by=c("concept_id"="measurement_concept_id")) %>% 
    compute_new(indexes=c("person_id"))
  
  if (test_only == TRUE) {
    rsv_tests %>% 
      filter(value_as_concept_id == 9191) %>% 
      # group_by(person_id) %>% 
      # slice_min(measurement_date, with_ties = FALSE) %>% 
      # ungroup() %>% 
      select(person_id, rsv_test_date = measurement_date,
             rsv_lab_concept_id = measurement_concept_id,
             rsv_lab_value_concept_id = value_as_concept_id) %>% 
      compute_new(indexes=c("person_id")) %>% 
      mutate(rsv_evidence_date = rsv_test_date) %>% 
      select(person_id, rsv_evidence_date,
             rsv_lab_concept_id,
             rsv_lab_value_concept_id) %>% 
      group_by(person_id, rsv_evidence_date) %>% 
      filter(row_number()==1) %>% 
      ungroup() %>% 
      return()
  } else {
    
    # Should be a full join
    rsv_tests %>% 
      filter(value_as_concept_id == 9191) %>% 
      # group_by(person_id) %>% 
      # slice_min(measurement_date, with_ties = FALSE) %>% 
      # ungroup() %>% 
      select(person_id, rsv_test_date = measurement_date) %>% 
      compute_new(indexes=c("person_id")) %>% 
      full_join(rsv_conditions %>% 
                  # group_by(person_id) %>% 
                  # # slice_min(condition_start_date, with_ties = FALSE) %>% # Same here, instead of slicing the min, could just identify the min
                  # ungroup() %>% 
                  select(person_id, rsv_dx_date = condition_start_date) %>% 
                  compute_new(indexes=c("person_id")), by=c("person_id")) %>% 
      mutate(rsv_evidence_date = case_when(is.na(rsv_test_date) ~ rsv_dx_date,
                                           is.na(rsv_dx_date) ~ rsv_test_date,
                                           rsv_dx_date <= rsv_test_date ~ rsv_dx_date,
                                           rsv_test_date < rsv_dx_date ~ rsv_test_date)) %>% 
      mutate(lab_confirmed_rsv = case_when(is.na(rsv_test_date) ~ 0,
                                           is.na(rsv_dx_date) ~ 1,
                                           rsv_dx_date <= rsv_test_date ~ 1,
                                           rsv_test_date < rsv_dx_date ~ 1)) %>% 
      select(person_id, rsv_evidence_date, lab_confirmed_rsv) %>% 
      group_by(person_id, rsv_evidence_date, lab_confirmed_rsv) %>% 
      filter(row_number()==1) %>% 
      ungroup() %>% 
      return()
  }

}

