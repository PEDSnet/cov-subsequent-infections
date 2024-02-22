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
  
  ### Step 1 
  ## Generate the 3 comparison cohorts
  ## Special inclusion/exclusion criteria for overlap in cohorts
  
  ###### Anchor dates for conditions
  postacute_min_start_date = as.Date("2022-04-01")
  postacute_max_end_date = as.Date("2023-07-01")
  
  cohort_entry_start_date = as.Date("2022-03-01")
  cohort_entry_end_date = as.Date("2023-01-01")

  cohort_1_label = "_rsv_study_cohort_"
  
  # odr <- cdm_tbl("observation_derivation_recover") %>% # Actually, need observations older than this, in order to get pre-index criteria
  #   filter(observation_date >= ce_start_date,
  #          observation_date < postacute_max_end_date)
  
  rslt$pats_under_age_6 <- get_patients_under_age_limit(cdm_tbl("person"),
                                                        cdm_tbl("visit_occurrence"), 
                                                        ce_start_date, 
                                                        ce_end_date, age_limit_years_rough = 6) %>% 
    compute_new(indexes=c("person_id"))
  
  ## TODO add to attrition step
  
  rslt$triple_cohort <- rslt$pats_under_age_6 %>% 
    build_comparison_cohorts_rsv_study(odr_tbl = odr, 
                                       ce_start_date = ce_start_date, 
                                       ce_end_date = ce_end_date) %>% 
    compute_new(indexes=c("person_id")) %>% 
    mutate(age_at_ce_under_limit = ifelse(age_years_on_ce_date < 5, 1, 0)) %>% 
    filter(age_at_ce_under_limit == 1)
  
  ## TODO add to attrition step
  
  rslt$cohort_inclusion_flags <- rslt$triple_cohort %>% 
    apply_inclusion_flags_rsv_study(ce_start_date = ce_start_date,
                                    ce_end_date = ce_end_date)
  
  ## TODO add to attrition step (study_eligible == 1)
  ## Attrition table also broken down by category / one table per category / add a column for the step type
  ## TODO output the table that's filtered to study_eligible but still has index based exclusion reasons
  
  rslt$index_infection_exclusions <- rslt$cohort_inclusion_flags %>% 
    filter(exclude == 0)
  
  ## TODO attrition step, based on detailed exclusion reason?
    
  ## Final cohort
  rslt$final_cohort <- rslt$cohort_inclusion_flags %>% 
    filter(study_eligible == 1) %>% 
    filter(exclude == 0)
  
  ## Step 2
  ## Generate outcomes
  rslt$final_cohort %>% 
    generate_outcome_dataset(rsv....) # TODO
  
  ## Step 3
  ## Generate covariates
  
  ## Step 4
  ## Data transformations
  
  
  ## Output final datasets at end, see next driver for modeling and/or reporting analytics!
  
  
  # By convention, accumulate execution results in a list rather than as
  # independent variables, in order to make returning the entire set easier
  rslt <- list()
  
  
  

  
  message('Done.')
  
  invisible(rslt)
  
}
