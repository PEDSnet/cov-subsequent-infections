# Vector of additional packages to load before executing the request
config_append('extra_packages', c())

#' Execute the request
#'
#' This function presumes the environment has been set up, and executes the
#' steps of the request.
#'
#' In addition to performing queries and analyses, the execution path in this
#' function should include periodic progress messages to the user, and logging
#' of intermediate totals and timing data through [append_sum()].
#'
#' @return The return value is dependent on the content of the request, but is
#'   typically a structure pointing to some or all of the retrieved data or
#'   analysis results.  The value is not used by the framework itself.
#' @md
.run  <- function() {

  setup_pkgs() # Load runtime packages as specified above

  message('Starting execution with framework version ',
          config('framework_version'))

  # Set up the step log with as many attrition columns as you need.
  # For example, this call sets up the log with a `persons` count that will be
  # required at each step.
  init_sum(cohort = 'Start', persons = 0)

  # By convention, accumulate execution results in a list rather than as
  # independent variables, in order to make returning the entire set easier
  rslt <- list()
  
  ###### Create cohort here:
  ce_start_date = as.Date("2021-04-01")
  ce_end_date = as.Date("2023-10-01")
  # cohort_1_label = "_06210223_"
  cohort_1_label = "_04211023_"
  odr <- cdm_tbl("observation_derivation_recover") %>% 
    filter(observation_date < ce_end_date, observation_date >= ce_start_date)
  
  rslt$pats_under_age_6 <- generate_output_unfiltered_cohort(
    ce_start_date = ce_start_date,
    ce_end_date = ce_end_date,
    cohort_label = cohort_1_label,
    odr = odr
  )
  
  # counts the number of persons for this part of the query
  append_sum(cohort = 'Patients in study period under age 6',
             persons = distinct_ct(rslt$pats_under_age_6))
  
  ## Retain and output demographic information for the unfiltered cohort,
  # in case we want to see if there are already any confounders in the dataset in terms of
  # exclusions we've made for utilization, age, etc. 
  unfiltered_cohort_demo <- generate_output_cohort_demo(
    ce_start_date = ce_start_date,
    ce_end_date = ce_end_date,
    cohort_tbl_name = paste0("cohort", cohort_1_label, "unfiltered"),
    cohort_output_tbl_name = paste0("cohort", cohort_1_label, "unfiltered_demo")
  )
  
  ## Next, apply washout logic to the cohort to obtain final sample size, 
  ## but one at a time with some attrition steps in between:
  cohort_with_washouts <-
    results_tbl(paste0("cohort", cohort_1_label, "unfiltered_demo")) %>% 
    mutate(sub_cohort = cohort_assignment) %>% 
    apply_washout_logic(odr = odr) %>% 
    compute_new()
  
  rslt$has_covid_or_flu <- cohort_with_washouts %>% 
    filter(cohort_assignment != "No infection")
  
  rslt$under_6_at_ce <- cohort_with_washouts %>% 
    filter(cohort_assignment != "No infection") %>% 
    filter(age_at_ce_under_limit == 1)
  
  rslt$has_visit_n_months_prior <- cohort_with_washouts %>% 
    filter(cohort_assignment != "No infection") %>% 
    filter(age_at_ce_under_limit == 1) %>% 
    filter(visit_within_n_years == 1)
  
  rslt$has_visit_in_follow_up <- cohort_with_washouts %>% 
    filter(cohort_assignment != "No infection") %>% 
    filter(age_at_ce_under_limit == 1) %>% 
    filter(visit_within_n_years == 1) %>% 
    filter(has_visit_in_followup_period == 1)
  
  rslt$has_no_other_washout <- cohort_with_washouts %>% 
    filter(cohort_assignment != "No infection") %>% 
    filter(age_at_ce_under_limit == 1) %>% 
    filter(visit_within_n_years == 1) %>% 
    filter(has_visit_in_followup_period == 1) %>% 
    filter(washout == 0)
  
  ## TODO determine more attrition steps as needed
  # rslt$eligible_with_influenza <- cohort_with_washouts %>% 
  #   filter(study_eligible == 1) %>% 
  #   filter(washout == 0) %>% 
  #   filter(cohort_assignment == "Influenza")
  # 
  # rslt$eligible_with_covid <- cohort_with_washouts %>% 
  #   filter(study_eligible == 1) %>% 
  #   filter(washout == 0) %>% 
  #   filter(cohort_assignment == "Covid")
    
  
  ## Attrition steps
  append_sum(cohort = 'Patients in study period under age 6 with evidence of SARS-CoV-2 or Influenza',
             persons = distinct_ct(rslt$has_covid_or_flu))
  
  append_sum(cohort = 'Patients in study period under age 6 at cohort entry date',
             persons = distinct_ct(rslt$under_6_at_ce))
  
  append_sum(cohort = 'Has a visit within 18 months of CE date',
             persons = distinct_ct(rslt$has_visit_n_months_prior))
  
  append_sum(cohort = 'Has a visit during the follow up period',
             persons = distinct_ct(rslt$has_visit_in_follow_up))
  
  append_sum(cohort = 'Does not fail other washout criteria',
             persons = distinct_ct(rslt$has_no_other_washout))
  
  
  ## TODO there may be other attrition steps we could add
  output_sum(name = "cohort_attrition")
  
  cohort_with_washouts %>% 
    filter(study_eligible == 1) %>% 
    filter(washout == 0) %>% 
    filter(cohort_assignment != "No infection") %>% 
    select(person_id, ce_date, sub_cohort = cohort_assignment,
           mock_medical_complexity, age_cat_at_ce, birth_date,
           sex_cat, cohort_entry_week, visit_type_admin_other_telemed,
           visit_type_outpatient, race_eth_cat, visit_type_inpatient_or_ed) %>% 
    output_tbl("base_cohort", indexes=c("person_id"))
  
  
  ### Now use the cohort to generate the analytic dataset with other variables
  
  

  message('Done.')

  invisible(rslt)

}
