#### Creating a cohort
### Inputs: cohort entry date minimum
### cohort entry date maximum
### cohort can later be edited down based on cohort entry date so better cast a wider net first

ce_start_date = as.Date("2021-06-01")
ce_end_date = as.Date("2023-02-01")
cohort_1_label = "_06210223_"
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

eligible_pats %>% 
  summarise_as_fraction(group = c(cohort_assignment), distinct_key = person_id)


