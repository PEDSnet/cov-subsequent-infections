#### Creating a cohort
### Inputs: cohort entry date minimum
### cohort entry date maximum
### cohort can later be edited down based on cohort entry date so better cast a wider net first

ce_start_date = as.Date("2021-04-01")
ce_end_date = as.Date("2023-10-01")
# cohort_1_label = "_06210223_"
cohort_1_label = "_04211023_"
odr <- cdm_tbl("observation_derivation_recover") %>% 
  filter(observation_date < ce_end_date, observation_date >= ce_start_date)

## Need to update and document cohort design
eligible_pats <-
  cdm_tbl("person") %>% 
  get_patients_under_age_limit(cdm_tbl("visit_occurrence"), ce_start_date, ce_end_date, age_limit_years_rough = 6) %>% 
  identify_cohort_group(odr_tbl = odr, ce_start_date = ce_start_date, ce_end_date = ce_end_date) %>% 
  compute_new(indexes=c("person_id")) %>% 
  mutate(age_at_ce_under_limit = ifelse(age_years_on_ce_date < 5, 1, 0)) %>% 
  mutate(index_date = ce_date) %>% 
  flag_visits_prior(cdm_tbl("visit_occurrence"), n_years = 1.5) %>% 
  flag_visit_follow_up(cdm_tbl("visit_occurrence")) %>% 
  select(-index_date) %>% 
  compute_new(indexes=c("person_id")) %>% 
  flag_study_eligiblity(ce_start_date, ce_end_date) %>% 
  compute_new(indexes=c("person_id"))

eligible_pats %>% 
  output_tbl(paste0("ce", cohort_1_label, "unfiltered"), indexes=c("person_id"))

results_tbl(paste0("ce", cohort_1_label, "unfiltered")) %>% 
  summarise_as_fraction(group = c(cohort_assignment), distinct_key = person_id)

eligible_pats_saved <- results_tbl(paste0("ce", cohort_1_label, "unfiltered"))
ce_start_date = as.Date("2021-04-01")
ce_end_date = as.Date("2023-10-01")
cohort_1_label = "_04211023_"

eligible_pats_demo <-
  eligible_pats_saved %>% 
  left_join(cdm_tbl("person") %>% 
              select(person_id, gender_concept_id,
                     race_concept_id, ethnicity_concept_id), by="person_id") %>% 
  mutate(mock_medical_complexity = case_when(person_id %% 17 == 1 ~ "complex chronic",
                                             person_id %% 8==1 ~ "chronic, non-complex",
                                             TRUE ~ "non complex or chronic")) %>% 
  mutate(sex_cat = case_when(gender_concept_id == 8507L ~ 'Male',
                             gender_concept_id == 8532L ~ 'Female',
                             TRUE ~ 'Other/unknown'),
         race_eth_cat = case_when(ethnicity_concept_id == 38003563L ~ 'Hispanic',
                                  race_concept_id == 8516L ~ 'NH_Black/AA',
                                  race_concept_id %in% c(8515L, 8557L) ~
                                    'NH_Asian/PI',
                                  #                               race_concept_id == 8657L ~ 'Native American',
                                  race_concept_id == 8527L ~ 'NH_White',
                                  race_concept_id == 44814659L ~ 'NH_Multiple',
                                  TRUE ~ 'Other/Unknown'),
         age_cat_at_ce = case_when(age_years_on_ce_date < 0.5 ~ "<6 months",
                                   age_years_on_ce_date < 1 ~ "6m-1y old",
                                   age_years_on_ce_date < 2 ~ "1-2y old",
                                   age_years_on_ce_date < 3 ~ "2-3y old",
                                   age_years_on_ce_date < 4 ~ "3-4y old",
                                   TRUE ~ "4-5y old")) %>% 
  mutate(cohort_entry_week = floor_date(ce_date, unit="week")) %>% 
  select(-gender_concept_id, -race_concept_id, -ethnicity_concept_id) %>% 
  flag_utilization_level(cdm_tbl("visit_occurrence"), lookback_days = 365) %>% 
  # left_join(results_tbl(paste0("cohort", cohort_1_label, "pmca")), by="person_id") %>% 
  compute_new(indexes=c("person_id"))

eligible_pats_demo %>% 
  output_tbl(paste0("ce", cohort_1_label, "unfiltered_demo"), indexes=c("person_id"))



