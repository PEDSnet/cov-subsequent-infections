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
  
  cohort_1_label = "resp_study"
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
  rslt$pats_under_age_6 %>% output_tbl(paste0(cohort_1_label, "_pats_under_6"))
  
  rslt$pats_under_age_6 <- results_tbl(paste0(cohort_1_label, "_pats_under_6"))
  
  # sample_pats_under_6 <- rslt$pats_under_age_6 %>% 
  #   slice_sample(n=3000) %>% 
  #   compute_new(indexes=c("person_id"))
  # 
  # sample_pats_under_6 %>% output_tbl("pat_sample_03_05")
  
  rslt$triple_cohort <- results_tbl(paste0(cohort_1_label, "pats_under_6")) %>%
    build_comparison_cohorts_resp_study_NN(odr_tbl = cdm_tbl("observation_derivation_recover"), # Should probably keep track of the inclusion concept code
                                           ce_start_date = cohort_entry_start_date, 
                                           ce_end_date = cohort_entry_end_date,
                                           cohort_1_label) %>% 
    compute_new(indexes=c("person_id")) %>% 
    mutate(age_at_ce_under_limit = ifelse(age_years_on_ce_date < 5, 1, 0)) %>% 
    filter(age_at_ce_under_limit == 1)
  rslt$triple_cohort %>% output_tbl(paste0(cohort_1_label, "_triple_cohort"))
  
  ## TODO add to attrition step
  ## TODO can filter by lab_confirmed flag or not
  rslt$cohort_inclusion_flags <- results_tbl(paste0(cohort_1_label, "_triple_cohort")) %>% 
    apply_inclusion_flags_rsv_study(ce_start_date = cohort_entry_start_date,
                                    ce_end_date = cohort_entry_end_date)
  rslt$cohort_inclusion_flags %>% output_tbl(paste0(cohort_1_label, "_cohort_inclusion_flags"))
  
  
  ## Attrition table also broken down by category / one table per category / add a column for the step type
  ## TODO output the table that's filtered to study_eligible but still has index based exclusion reasons
  
  ## TODO attrition step, based on detailed exclusion reason?
  results_tbl(paste0(cohort_1_label, "_cohort_inclusion_flags")) %>% 
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
  
  ## Save table with exclusion reasons
  rslt$other_washout_reasons %>% 
    filter(study_eligible == 1) %>% 
    output_tbl(paste0(cohort_1_label, "_w_excl_reasons"), indexes=c("person_id"))
  
  ## Final cohort
  ## CORRECT EXCLUSION REASONS
  rslt$final_cohort <- rslt$other_washout_reasons %>% 
    mutate(exclude = case_when(exclude_reason == "Had COVID: post-acute, FLU: index" ~ 0,
                               exclude_reason == "Had COVID: index, FLU: post-acute" ~ 0,
                               TRUE ~ exclude)) %>% 
    filter(study_eligible == 1) %>% 
    filter(exclude != 1) %>% 
    filter(washout == 0)
  
  ## Save table with exclusion reasons
  rslt$final_cohort %>% 
    output_tbl(paste0(cohort_1_label, "_cohort"))
  
  rslt$final_cohort_demo <-
    rslt$final_cohort %>% 
    cohort_demo() %>% 
    get_insurance_class() 
  
  rslt$final_cohort_demo %>% 
    output_tbl(paste0(cohort_1_label, "_cohort"), indexes=c("person_id"))
  
  ## Step 2
  ## Generate outcomes
  results_tbl(paste0(cohort_1_label, "_cohort")) %>% 
    generate_output_outcome_lists(outcome = "respiratory",
                                  outcome_start_date = cohort_entry_start_date,
                                  outcome_end_date = postacute_max_end_date,
                                  output_tbl_name = paste0(cohort_1_label, "_resp_outcomes"))
  ## In addition, get covid outcomes as well, or record from pre-existing data to get a patient's first covid or reinfection outcome
  
  ## Generate outcomes (any infection)
  results_tbl(paste0(cohort_1_label, "_cohort"), results_tag = FALSE) %>% 
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
  cohort_with_util_info <- results_tbl(paste0(cohort_1_label, "_cohort")) %>% 
    # select(-visit_type_outpatient, -visit_type_admin_other_telemed, -visit_type_inpatient_or_ed) %>% # remove previous bad versions 
    flag_utilization_level(cdm_tbl("visit_occurrence"), 365, cohort_1_label)
  
  cohort_with_util_info %>% 
    select(person_id, visit_type_inpatient, visit_type_outpatient, visit_type_ed, visit_type_other) %>% 
    output_tbl(paste0(cohort_1_label, "_visit_quantile_cat"))
  
  ## Hospitalization
  cohort_with_hosp_info <- results_tbl(paste0(cohort_1_label, "_cohort")) %>% 
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
    output_tbl(paste0(cohort_1_label, "_severe_ventilator"))
  
  ### PMCA 
  ##############################################################################
  message('Steps 1 and 2: Calculate PMCA lookup table and summary table for all sites in cohort')
  
  #' The specified cohort should have person_id and observation_date; a 3-year lookback from
  #' the observation date is applied
  pmca_lookup <- produce_pmca_lookup(cohort= results_tbl(paste0(cohort_1_label, "_cohort")) %>% 
                                       select(-observation_date) %>% 
                                       mutate(observation_date = ce_date))
  pmca_lookup %>% 
    output_tbl(paste0(cohort_1_label, '_pmca_lookup'))
  
  pmca_summary <- compute_pmca_summary(pmca_lookup_tbl = pmca_lookup)
  pmca_summary %>% output_tbl(paste0(cohort_1_label, '_pmca_summary'))
  
  ##############################################################################
  message('Step 3: Apply most conservative PMCA algorithm and get PMCA flags')
  
  pmca_category <- compute_pmca_cats_cons(pmca_summary_tbl = pmca_summary,
                                          cohort_tbl = results_tbl(paste0(cohort_1_label, "_cohort"))  %>% 
                                            inner_join(cdm_tbl("person") %>% select(person_id, site), by="person_id"))
  
  output_tbl(pmca_category,
             paste0(cohort_1_label, '_pmca_category'),
             indexes=c('person_id'))
  
  ##############################################################################
  message('Steps 4 and 5: Compute PMCA complex chronic/progressive/malignant and body system flags')
  
  #' Compute flags, including:
  #' progressive_ct
  #' malignancy_ct
  #' complex_chronic, chronic, non_complex_chronic
  #' n_body_systems
  pmca_all_flags_tbl <- pmca_all_flags(pmca_summary_tbl = pmca_summary,
                                       cohort_tbl = results_tbl(paste0(cohort_1_label, "_cohort")) ,
                                       pmca_category = pmca_category)
  pmca_all_flags_tbl %>% 
    output_tbl(paste0(cohort_1_label, '_pmca_all_flags'),
               indexes=c('person_id'))
  #' Compute flags for any indication of each PMCA body system in the prior 3 years
  pmca_bs_flags <- flag_pmca_bs(cohort_tbl= results_tbl(paste0(cohort_1_label, "_cohort")) ,
                                pmca_summary=pmca_summary)
  pmca_bs_flags %>% 
    output_tbl(paste0(cohort_1_label, '_pmca_bs_flags_full'),
               indexes=c('person_id'))
  ##############################################################################
  message('Define PMCA flags in your cohort')
  #' Note: the complex_chronic, non_complex_chronic, and chronic flags are oddly defined:
  #' complex_chronic: flags patients who are identified as complex AND chronic patients
  #' chronic: flags patients who are identified as chronic patients (but not complex)
  #' non_complex_chronic: flags patients who are neither "complex chronic" nor "chronic" patients
  #' (patients who do not have any indication of chronic or complex chronic disease)
  
  cohort_pmca_flags <- results_tbl(paste0(cohort_1_label, "_cohort")) %>%
    select('person_id') %>% 
    left_join(results_tbl(paste0(cohort_1_label, '_pmca_all_flags')) %>%
                select(-site), by='person_id') %>%
    left_join(results_tbl(paste0(cohort_1_label, '_pmca_bs_flags_full')), by='person_id') %>%
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
  
  
  ## Step 4
  ## Data transformations
  
  ## Output final analytic dataset at end, see next driver for modeling and/or reporting analytics!
  ## TODO have this be a running / evolving piece of code that generates the "final analytic dataset" that is ready to be plugged in to different models
  results_tbl(paste0(cohort_1_label, "_cohort")) %>% 
    #select(-visit_type_outpatient, -visit_type_admin_other_telemed, -visit_type_inpatient_or_ed) %>% # remove previous bad versions 
    generate_analytic_dataset_resp(respiratory_outcomes = results_tbl(paste0(cohort_1_label, "_resp_outcomes")),
                                   any_infection_outcomes = results_tbl(paste0(cohort_1_label, "_any_infection_outcomes")),
                                   util_quantiles_tbl = results_tbl(paste0(cohort_1_label, "_visit_quantile_cat")),
                                   hospitalization_tbl = results_tbl(paste0(cohort_1_label, "_hospitalization")),
                                   ventilator_use_tbl = results_tbl(paste0(cohort_1_label, "_severe_ventilator")),
                                   pmca_tbl = results_tbl(paste0(cohort_1_label, "_pmca_info")),
                                   output_tbl_name = paste0(cohort_1_label, "_analytic_dataset"))
  
  ## Step 5, get high risk flags
  analytic_dataset <-
    results_tbl(paste0(cohort_1_label, "_analytic_dataset")) %>% compute_new()
  
  highrisk_tbl <- analytic_dataset %>% select(person_id, ce_date) %>% 
    left_join(cdm_tbl("condition_occurrence") %>% 
                select(condition_source_concept_id, person_id, condition_start_date),# %>% compute_new(indexes=c("person_id", "condition_source_concept_id")), 
              by="person_id") %>% 
    # filter(ce_date >= condition_start_date, ce_date <= condition_start_date + days(1096))%>%# 3 years look-back for highrisk
    distinct(person_id, condition_source_concept_id, .keep_all = TRUE)%>%
    compute_new(indexes=c("person_id", "condition_source_concept_id"))
  
  dx_highrisk = load_codeset("dx_highrisk")
  highrisk_tbl <- highrisk_tbl %>% inner_join(dx_highrisk %>% select(concept_id) %>% mutate(highrisk_flag = TRUE),
                                              by = c("condition_source_concept_id" = "concept_id")) %>%
    select(-condition_source_concept_id) %>%
    distinct()%>% compute_new()
  highrisk_tbl %>% output_tbl(paste0(cohort_1_label, "_highrisk_dx_1"))
  
  highrisk_tbl <- results_tbl(paste0(cohort_1_label, "_highrisk_dx_1"))
  analytic_dataset_highrisk <- analytic_dataset %>% 
    left_join(results_tbl(paste0(cohort_1_label, "_highrisk_dx_1")) %>% 
                # filter(condition_start_date <= ce_date, 
                #         ce_date - condition_start_date >= 1096)%>%
                group_by(person_id) %>% 
                slice_max(condition_start_date, n = 1, with_ties = FALSE) %>% ungroup() %>%
                select(-ce_date), by = "person_id") %>%
    mutate(highrisk_flag = if_else(is.na(highrisk_flag), FALSE, TRUE))
  analytic_dataset_highrisk %>% output_tbl(paste0(cohort_1_label, "_analytic_highrisk_1"))
  
  ## Now get the 4x4 table that Suchitra wants
  analytic_dataset_highrisk <- results_tbl(paste0(cohort_1_label, "_analytic_highrisk_1")) 
  analytic_dataset_highrisk %>% distinct(ce_event, lab_confirmed, sub_cohort)
  
  dx_only <- analytic_dataset_highrisk %>% filter(sub_cohort == "Influenza")%>%
    anti_join(results_tbl(paste0(cohort_1_label, "_flu_lab_complete")) %>%
                distinct(person_id), by = "person_id") %>%
    distinct(person_id) %>%
    mutate(dx_only = TRUE) %>% compute_new()
  # dx_only %>% distinct_ct()
  
  lab_only <- analytic_dataset_highrisk %>% filter(sub_cohort == "Influenza")%>%
    anti_join(results_tbl(paste0(cohort_1_label, "_flu_dx"))%>%
                distinct(person_id), by = "person_id") %>%
    distinct(person_id) %>%
    mutate(lab_only = TRUE) %>% compute_new()
  # lab_only %>% distinct_ct()
  
  dx_lab <- analytic_dataset_highrisk %>% filter(sub_cohort == "Influenza")%>%
    anti_join(lab_only, by = "person_id") %>%
    anti_join(dx_only, by = "person_id") %>%
    distinct(person_id) %>%
    mutate(dx_lab = TRUE) %>% compute_new()
  # dx_lab %>% distinct_ct()
  
  analytic_dataset_highrisk %>% left_join(lab_only, by = "person_id") %>%
    left_join(dx_only, by = "person_id") %>%
    left_join(dx_lab, by = "person_id") %>%
    mutate(flu_type = case_when(dx_lab ~ "dx_lab", dx_only ~ "dx_only", lab_only ~ "lab_only", TRUE ~NA))%>%
    # rename(ce_date = ce_date.x) %>%
    select(-dx_only, -lab_only, -dx_lab) %>%
    output_tbl(paste0(cohort_1_label, "_analytic_final_1"))
  
  message("output attritrion table")
  
  message("sites with < 20% confirmed influenza test are excluded")
  site_excl <- c("monte", "nationwide", "columbia", "emory", 
                 "intermountain", "ochsner", "temple", "uth", 
                 "osu", "utsw", "wakeforest")
  
  output_attrition_table_site_excl(cohort_1_label = "resp_study", 
                         cohort_entry_start_date = as.Date("2021-06-01"),
                         cohort_entry_end_date = as.Date("2023-01-01"),
                         site_excl)
  
  message('Done.')
  
  invisible(rslt)
  
}
