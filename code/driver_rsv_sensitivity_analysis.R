## Specific driver for the portion of the study limited to the RSV outcome / comparitive risk study

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
  
  init_sum(cohort = 'Start', persons = 0)
  
  rslt = list()
  
  ### Step 1 
  ## Generate the 3 comparison cohorts
  ## Special inclusion/exclusion criteria for overlap in cohorts
  
  ###### Anchor dates for conditions
  postacute_min_start_date = as.Date("2022-04-01")
  postacute_max_end_date = as.Date("2023-07-01")
  
  cohort_entry_start_date = as.Date("2022-03-01")
  cohort_entry_end_date = as.Date("2023-01-01")
  
  cohort_1_label = "_rsv_study_cohort_"
  cohort_1_label = "rsv_sensitivity_"
  
  # odr <- cdm_tbl("observation_derivation_recover") %>% # Actually, need observations older than this, in order to get pre-index criteria
  #   filter(observation_date >= ce_start_date,
  #          observation_date < postacute_max_end_date)
  
  ## Easier way to do this is just select patients who have a ce_end_date - date of birth < 6 
  rslt$pats_under_age_6 <- pats_with_visit_under_age(cdm_tbl("person"),
                                                     cdm_tbl("visit_occurrence"), 
                                                     cohort_entry_start_date, 
                                                     cohort_entry_end_date, age_limit_years_rough = 6) %>% 
    compute_new(indexes=c("person_id"))
  
  append_sum(cohort = 'Patients under 6 years old with at least 1 visit during CE period',
             persons = distinct_ct(rslt$pats_under_age_6))
  
  sample_pats_under_6 <- rslt$pats_under_age_6 %>% 
    slice_sample(n=500) %>% 
    compute_new(indexes=c("person_id"))
  
  sample_pats_under_6 %>% output_tbl("pat_sample_02_22")
  
  rslt$triple_cohort_all <- rslt$pats_under_age_6 %>% 
    build_comparison_cohorts_rsv_study_sensitivity(odr_tbl = cdm_tbl("observation_derivation_recover"), # Should probably keep track of the inclusion concept code
                                       ce_start_date = cohort_entry_start_date, 
                                       ce_end_date = cohort_entry_end_date) %>% 
    compute_new(indexes=c("person_id")) 
  
  rslt$ triple_cohort <- rslt$triple_cohort_all %>% 
    mutate(age_at_ce_under_limit = ifelse(age_years_on_ce_date < 5, 1, 0)) %>% 
    filter(age_at_ce_under_limit == 1)
  
  ## TODO add to attrition step
  append_sum(cohort = 'Had at least 1 of COVID, influenza, and/or other respiratory infection event during CE period',
             persons = distinct_ct(rslt$triple_cohort))
  
  append_sum(cohort = 'Had COVID during CE period',
             persons = distinct_ct(rslt$triple_cohort %>% 
                                     filter(!is.na(covid_date),
                                            covid_date >= cohort_entry_start_date,
                                            covid_date < cohort_entry_end_date)))
  
  append_sum(cohort = 'Had influenza during CE period',
             persons = distinct_ct(rslt$triple_cohort %>% 
                                     filter(!is.na(flu_date),
                                            flu_date >= cohort_entry_start_date,
                                            flu_date < cohort_entry_end_date)))
  
  append_sum(cohort = 'Had other respiratory infection during CE period',
             persons = distinct_ct(rslt$triple_cohort %>% 
                                     filter(!is.na(resp_date),
                                            resp_date >= cohort_entry_start_date,
                                            resp_date < cohort_entry_end_date)))
  
  
  rslt$cohort_inclusion_flags <- rslt$triple_cohort %>% 
    apply_inclusion_flags_rsv_study(ce_start_date = cohort_entry_start_date,
                                    ce_end_date = cohort_entry_end_date)
  
  ## TODO add to attrition step (study_eligible == 1)
  append_sum(cohort = 'Meets utilization inclusion criteria',
             persons = distinct_ct(rslt$cohort_inclusion_flags %>% 
                                     filter(study_eligible == 1)))
  
  append_sum(cohort = 'Not excluded due to overlapping index infections',
             persons = distinct_ct(rslt$cohort_inclusion_flags %>% 
                                     filter(study_eligible == 1,
                                            exclude == 0)))
  
  ## Attrition table also broken down by category / one table per category / add a column for the step type
  ## TODO output the table that's filtered to study_eligible but still has index based exclusion reasons
  
  ## TODO attrition step, based on detailed exclusion reason?
  rslt$cohort_inclusion_flags %>% 
    filter(study_eligible == 1) %>% 
    group_by(sub_cohort, exclude_reason) %>% 
    summarise(n=n_distinct(person_id)) %>% 
    output_tbl(paste0(cohort_1_label, "_exclusion_reasons"))
  
  rslt$other_washout_reasons <-
    rslt$cohort_inclusion_flags %>% 
    apply_washout_logic(odr = cdm_tbl("observation_derivation_recover")) %>% 
    compute_new()
  
  rslt$other_washout_reasons %>% 
    filter(study_eligible == 1) %>% 
    filter(exclude == 0) %>% 
    group_by(sub_cohort, washout_reason) %>% 
    summarise(n=n_distinct(person_id)) %>% 
    output_tbl(paste0(cohort_1_label, "_washout_reasons"))
  
  ## Final cohort
  rslt$final_cohort <- rslt$other_washout_reasons %>% 
    filter(study_eligible == 1) %>% 
    filter(exclude == 0) %>% 
    filter(washout == 0)
  
  ## Save table with exclusion reasons
  rslt$other_washout_reasons %>% 
    filter(study_eligible == 1) %>% 
    output_tbl(paste0(cohort_1_label, "_w_exclusion_reasons"), indexes=c("person_id"))
  
  append_sum(cohort = 'Final COVID sub-cohort',
             persons = distinct_ct(rslt$final_cohort %>% 
                                     filter(sub_cohort=="COVID")))
  
  append_sum(cohort = 'Final Influenza sub-cohort',
             persons = distinct_ct(rslt$final_cohort %>% 
                                     filter(sub_cohort=="Influenza")))
  
  append_sum(cohort = 'Final Other respiratory sub-cohort',
             persons = distinct_ct(rslt$final_cohort %>% 
                                     filter(sub_cohort=="Respiratory")))
  
  output_sum(name = paste0(cohort_1_label, "_attrition"))
  
  rslt$final_cohort_demo <-
    rslt$final_cohort %>% 
    cohort_demo() %>% 
    get_insurance_class() 
  
  rslt$final_cohort_demo %>% 
    output_tbl(paste0(cohort_1_label, "_cohort_demo"), indexes=c("person_id"))
  
  ## Step 2
  ## Generate outcomes
  results_tbl(paste0(cohort_1_label, "_cohort_demo")) %>% 
    generate_output_outcome_lists(outcome = "rsv",
                                  outcome_start_date = as.Date("2022-04-01"),
                                  outcome_end_date = as.Date("2023-07-01"),
                                  output_tbl_name = paste0(cohort_1_label, "rsv_outcomes"),
                                  test_only = TRUE)
  
  ## FIRST THING TUESDAY: Run report using the results of the sensitivity analysis, see what the distributions show
  
  ## Step 3
  ## Generate additional covariates
  
  ## Step 4
  ## Data transformations
  
  
  ## Output final datasets at end, see next driver for modeling and/or reporting analytics!
  
  
  
  message('Done.')
  
  invisible(rslt)
  
}
