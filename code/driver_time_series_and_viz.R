
postacute_min_start_date = as.Date("2022-04-01")
postacute_max_end_date = as.Date("2023-01-01")

cohort_entry_start_date = as.Date("2022-03-01")
cohort_entry_end_date = as.Date("2022-07-01")

cohort_1_label = "rsv_reverse_analysis"

pats_under_age_6 <- cdm_tbl("person") %>% 
  filter(birth_date > as.Date("2017-07-01")) %>% 
  select(person_id, site) %>% 
  compute_new(indexes=c("site", "person_id"))

all_rsv_cases_in_outcome_period <- pats_under_age_6 %>% 
  flag_rsv_outcome(outcome_start_date = as.Date("2022-04-01"),
                   outcome_end_date = as.Date("2023-01-01")
  )

# Then: What fraction of them had covid and what portion of them had flu in the 30-180 days before, or something else / none
all_rsv_cases_in_outcome_period %>% 
  group_by(person_id) %>% 
  slice_min(rsv_evidence_date, with_ties = FALSE) %>% 
  ungroup() %>% 
  mutate(lc=as.character(lab_confirmed_rsv)) %>%
  ggplot() + 
  geom_density(aes(x=rsv_evidence_date, fill=lc), position="stack")

all_rsv_cases_in_outcome_period %>% 
  output_tbl("all_rsv_cases")

all_rsv_with_covid <- results_tbl("all_rsv_cases") %>% 
  group_by(person_id) %>% 
  slice_min(rsv_evidence_date, with_ties = FALSE) %>% 
  ungroup() %>% 
  mutate(min_index_date = as.Date(rsv_evidence_date - days(180))) %>% 
  mutate(max_index_date = as.Date(rsv_evidence_date - days(30))) %>% 
  flag_covid_positives(cdm_tbl("observation_derivation_recover")) %>% 
  compute_new(indexes=c("person_id"))

all_rsv_with_covid %>% 
  mutate(had_covid_in_window = case_when(
    covid_event_positive==1 & earliest_cov_event_date <= max_index_date & earliest_cov_event_date >= min_index_date ~ 1,
    TRUE ~ 0
  )) %>% 
  mutate(had_covid_in_window = as.character(had_covid_in_window)) %>% 
  ggplot() +
  geom_density(aes(x=rsv_evidence_date, fill = had_covid_in_window), position="stack")

all_rsv_with_flu <- results_tbl("all_rsv_cases") %>% 
  group_by(person_id) %>% 
  slice_min(rsv_evidence_date, with_ties = FALSE) %>% 
  ungroup() %>% 
  mutate(min_index_date = as.Date(rsv_evidence_date - days(180))) %>% 
  mutate(max_index_date = as.Date(rsv_evidence_date - days(30))) %>% 
  get_flu_evidence(min_date = cohort_entry_start_date, max_date = cohort_entry_end_date) %>% 
  compute_new(indexes=c("person_id"))

full_disease_panel <- all_rsv_with_covid %>% 
  select(person_id, rsv_evidence_date, min_index_date, max_index_date, earliest_cov_event_date, covid_event_positive) %>% 
  left_join(all_rsv_with_flu, by="person_id") %>% 
  mutate(had_covid_in_window = case_when(
    covid_event_positive==1 & earliest_cov_event_date <= max_index_date & earliest_cov_event_date >= min_index_date ~ 1,
    TRUE ~ 0
  )) %>% 
  mutate(had_flu_in_window = case_when(
    earliest_flu_evidence <= max_index_date & earliest_flu_evidence >= min_index_date ~ 1,
    TRUE ~ 0
  )) %>% 
  mutate(pre_acute_disease = case_when(
    had_covid_in_window==1 & had_flu_in_window==1 ~ "Both",
    had_covid_in_window==1 ~ "SARS-CoV-2",
    had_flu_in_window==1 ~ "Influenza",
    TRUE ~ "None"
  )) %>% 
  compute_new()

full_disease_panel %>% 
  output_tbl("descriptive_rsv_surge_data")
  



full_disease_panel %>% 
  ggplot() +
  geom_density(aes(x=rsv_evidence_date, fill=pre_acute_disease), position="stack")

table1(~ pre_acute_disease, data = full_disease_panel)

## Make figure here for the table of the full landscape of infections
# Use the existing cohort
# Use the new RSV disease panel
# Join
# Make a figure

### NHAT HERE IS THE CODE FOR THE PRETTY FIGURE THAT SUCHITRA WANTS!

results_tbl("resp_study_cohort_analytic_dataset") %>% 
  select(person_id, infection_date = ce_date, disease = sub_cohort) %>% 
  union(results_tbl("descriptive_rsv_surge_data") %>% 
          select(person_id, infection_date = rsv_evidence_date) %>% 
          mutate(disease = "rsv")) %>% 
  mutate(infection_week = floor_date(infection_date, unit="week")) %>% 
  group_by(infection_week, disease) %>% 
  summarise(n=n_distinct(person_id)) %>% 
  ungroup() %>% 
  ggplot() +
  geom_bar(aes(x=infection_week, y=n, fill=disease), position="stack", stat="identity") +
  scale_fill_brewer(palette = "Dark2") +
  theme_bw()









