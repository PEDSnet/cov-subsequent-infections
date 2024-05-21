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
  
  ## Can we use the same cohorts for the next analysis, just change the preclusion criteria and add additional outcomes?
  ###### Anchor dates for conditions
  postacute_min_start_date = as.Date("2021-07-01")
  postacute_max_end_date = as.Date("2023-03-01")
  
  cohort_entry_start_date = as.Date("2021-06-01")
  cohort_entry_end_date = as.Date("2023-02-01")
  
  cohort_1_label = "resp_study_cohort"
  # cohort_1_label = "sample_pats"
  
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
  
  # sample_pats_under_6 <- rslt$pats_under_age_6 %>% 
  #   slice_sample(n=3000) %>% 
  #   compute_new(indexes=c("person_id"))
  # 
  # sample_pats_under_6 %>% output_tbl("pat_sample_03_05")
  
  rslt$triple_cohort <- rslt$pats_under_age_6 %>% 
    build_comparison_cohorts_resp_study(odr_tbl = cdm_tbl("observation_derivation_recover"), # Should probably keep track of the inclusion concept code
                                       ce_start_date = cohort_entry_start_date, 
                                       ce_end_date = cohort_entry_end_date) %>% 
    compute_new(indexes=c("person_id")) %>% 
    mutate(age_at_ce_under_limit = ifelse(age_years_on_ce_date < 5, 1, 0)) %>% 
    filter(age_at_ce_under_limit == 1)
  
  ## TODO add to attrition step
  ## TODO can filter by lab_confirmed flag or not
  append_sum(cohort = 'Had at least 1 of COVID, influenza infection event during CE period',
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

  
  rslt$cohort_inclusion_flags <- rslt$triple_cohort %>% 
    apply_inclusion_flags_rsv_study(ce_start_date = cohort_entry_start_date,
                                    ce_end_date = cohort_entry_end_date)
  
  ## TODO add to attrition step (study_eligible == 1)
  append_sum(cohort = 'Meets utilization inclusion criteria',
             persons = distinct_ct(rslt$cohort_inclusion_flags %>% 
                                     filter(study_eligible == 1)))
  
  
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
    filter(exclude != 1) %>% 
    group_by(sub_cohort, washout_reason) %>% 
    summarise(n=n_distinct(person_id)) %>% 
    output_tbl(paste0(cohort_1_label, "_washout_reasons"))
  
  ## Final cohort
  ## CORRECT EXCLUSION REASONS
  rslt$final_cohort <- rslt$other_washout_reasons %>% 
    mutate(exclude = case_when(exclude_reason == "Had COVID: post-acute, FLU: index" ~ 0,
                               exclude_reason == "Had COVID: index, FLU: post-acute" ~ 0,
                               TRUE ~ exclude)) %>% 
    filter(study_eligible == 1) %>% 
    filter(exclude != 1) %>% 
    filter(washout == 0)
  
  append_sum(cohort = 'Not excluded due to overlapping index infections',
             persons = distinct_ct(rslt$final_cohort %>% 
                                     filter(study_eligible == 1,
                                            exclude != 1)))
  
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
    generate_output_outcome_lists(outcome = "respiratory",
                                  outcome_start_date = cohort_entry_start_date,
                                  outcome_end_date = postacute_max_end_date,
                                  output_tbl_name = paste0(cohort_1_label, "_respiratory_outcomes"))
  ## In addition, get covid outcomes as well, or record from pre-existing data to get a patient's first covid or reinfection outcome
  
  ## Generate outcomes (any infection)
  results_tbl(paste0(cohort_1_label, "_cohort_demo")) %>% 
    generate_output_outcome_lists(outcome = "general",
                                  outcome_start_date = cohort_entry_start_date,
                                  outcome_end_date = postacute_max_end_date,
                                  output_tbl_name = paste0(cohort_1_label, "_any_infection_outcomes"))
  
  ## Step 3
  ## Generate additional covariates
  ## In person visit around index (aka severity)
  ## Body systems
  
  
  ### Utilization here
  ## TODO continue this on 3/26
  cohort_with_util_info <- results_tbl(paste0(cohort_1_label, "_cohort_demo")) %>% 
    # select(-visit_type_outpatient, -visit_type_admin_other_telemed, -visit_type_inpatient_or_ed) %>% # remove previous bad versions 
    flag_utilization_level(cdm_tbl("visit_occurrence"), 365, cohort_1_label)
  
  cohort_with_util_info %>% 
    select(person_id, visit_type_inpatient, visit_type_outpatient, visit_type_ed, visit_type_other) %>% 
    output_tbl(paste0(cohort_1_label, "_visit_quantile_cat"))
  
  ## Hospitalization
  cohort_with_hosp_info <- results_tbl(paste0(cohort_1_label, "_cohort_demo")) %>% 
    mutate(index_date = ce_date) %>% 
    flag_hospitalizations()
  
  cohort_with_hosp_info %>% 
    select(person_id, hospital_flag, hospitalization_visit_date, visit_span_criteria_rank ) %>% 
    output_tbl(paste0(cohort_1_label, "_hospitalization"))
  
  ## Ventilator use
  cohort_with_vent_info <- cohort_with_hosp_info %>% 
    # mutate(index_date = ce_date) %>% 
    get_severe_hospitalizations()
  
  cohort_with_vent_info %>% 
    select(person_id, earliest_proc_vent_date, hosp_with_vent ) %>% 
    output_tbl(paste0(cohort_1_label, "_severe_ventilator_use"))
  
  ### PMCA 
  ##############################################################################
  message('Steps 1 and 2: Calculate PMCA lookup table and summary table for all sites in cohort')
  
  #' The specified cohort should have person_id and observation_date; a 3-year lookback from
  #' the observation date is applied
  pmca_lookup <- produce_pmca_lookup(cohort= results_tbl(paste0(cohort_1_label, "_cohort_demo")) %>% 
                                       select(-observation_date) %>% 
                                       mutate(observation_date = ce_date))
  
  pmca_summary <- compute_pmca_summary(pmca_lookup_tbl = pmca_lookup)
  
  output_tbl(pmca_summary,
             paste0(cohort_1_label,'pmca_summary'),
             indexes=c('person_id'))
  
  ##############################################################################
  message('Step 3: Apply most conservative PMCA algorithm and get PMCA flags')
  
  pmca_category <- compute_pmca_cats_cons(pmca_summary_tbl = pmca_summary,
                                          cohort_tbl = results_tbl(paste0(cohort_1_label, "_cohort_demo"))  %>% 
                                            inner_join(cdm_tbl("person") %>% select(person_id, site), by="person_id"))
  
  # output_tbl(rslt$pmca_category,
  #            'pmca_category',
  #            indexes=c('person_id'))
  
  ##############################################################################
  message('Steps 4 and 5: Compute PMCA complex chronic/progressive/malignant and body system flags')
  
  #' Compute flags, including:
  #' progressive_ct
  #' malignancy_ct
  #' complex_chronic, chronic, non_complex_chronic
  #' n_body_systems
  pmca_all_flags <- pmca_all_flags(pmca_summary_tbl = pmca_summary,
                                   cohort_tbl = results_tbl(paste0(cohort_1_label, "_cohort_demo")) ,
                                   pmca_category = pmca_category)
  
  #' Compute flags for any indication of each PMCA body system in the prior 3 years
  pmca_bs_flags <- flag_pmca_bs(cohort_tbl= results_tbl(paste0(cohort_1_label, "_cohort_demo")) ,
                                pmca_summary=pmca_summary)
  
  ##############################################################################
  message('Define PMCA flags in your cohort')
  #' Note: the complex_chronic, non_complex_chronic, and chronic flags are oddly defined:
  #' complex_chronic: flags patients who are identified as complex AND chronic patients
  #' chronic: flags patients who are identified as chronic patients (but not complex)
  #' non_complex_chronic: flags patients who are neither "complex chronic" nor "chronic" patients
  #' (patients who do not have any indication of chronic or complex chronic disease)
  
  cohort_pmca_flags <- results_tbl(paste0(cohort_1_label, "_cohort_demo")) %>%
    select('person_id') %>% 
    left_join(pmca_all_flags %>%
                select(-site), by='person_id') %>%
    left_join(pmca_bs_flags, by='person_id') %>%
    mutate(complex_chronic_flag = case_when(complex_chronic == 1L ~ 'complex_chronic',
                                            non_complex_chronic == 1L ~ 'not_chronic_or_complex_chronic',
                                            chronic == 1L ~ 'chronic',
                                            TRUE ~ 'not_chronic_or_complex_chronic')) %>%
    mutate(progressive_flag = case_when(progressive_ct > 0L ~ 'yes',
                                        TRUE ~ 'no'),
           malignancy_flag = case_when(malignancy_ct > 0L ~ 'yes',
                                       TRUE ~ 'no')) %>%
    mutate(n_body_systems=case_when(is.na(n_body_systems) ~ 0L,
                                    TRUE ~ n_body_systems)) %>%
    mutate(nbs=case_when(
      n_body_systems>=5L~"05-17 body systems",
      n_body_systems>=3L~"03-04 body systems",
      n_body_systems==2L~"02 body systems",
      n_body_systems==1L~"01 body system",
      n_body_systems==0L~"No body systems")) %>% 
    compute_new()
  
  cohort_pmca_flags %>% 
    output_tbl(paste0(cohort_1_label, "_pmca_info"))
  
  # output_tbl(rslt$pmca_bs_flags,
  #            'pmca_bs_flags_full',
  #            indexes=list('person_id'))
  
  
  ## Step 4
  ## Data transformations
  
  ## Output final analytic dataset at end, see next driver for modeling and/or reporting analytics!
  ## TODO have this be a running / evolving piece of code that generates the "final analytic dataset" that is ready to be plugged in to different models
  results_tbl(paste0(cohort_1_label, "_cohort_demo")) %>% 
    # select(-visit_type_outpatient, -visit_type_admin_other_telemed, -visit_type_inpatient_or_ed) %>% # remove previous bad versions 
    generate_analytic_dataset_resp(respiratory_outcomes = results_tbl(paste0(cohort_1_label, "_respiratory_outcomes")),
                                   any_infection_outcomes = results_tbl(paste0(cohort_1_label, "_any_infection_outcomes")),
                              util_quantiles_tbl = results_tbl(paste0(cohort_1_label, "_visit_quantile_cat")),
                              hospitalization_tbl = results_tbl(paste0(cohort_1_label, "_hospitalization")),
                              ventilator_use_tbl = results_tbl(paste0(cohort_1_label, "_severe_ventilator_use")),
                              pmca_tbl = results_tbl(paste0(cohort_1_label, "_pmca_info")),
                              output_tbl_name = paste0(cohort_1_label, "_analytic_dataset"))
  
  message('Done.')
  
  invisible(rslt)
  
}
