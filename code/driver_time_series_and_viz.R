
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
# all_rsv_cases_in_outcome_period %>% 
#   group_by(person_id) %>% 
#   slice_min(rsv_evidence_date, with_ties = FALSE) %>% 
#   ungroup() %>% 
#   mutate(lc=as.character(lab_confirmed_rsv)) %>%
#   ggplot() + 
#   geom_density(aes(x=rsv_evidence_date, fill=lc), position="stack")

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

all_rsv_with_covid %>% output_tbl("all_rsv_w_covid")
# all_rsv_with_covid %>% 
#   mutate(had_covid_in_window = case_when(
#     covid_event_positive==1 & earliest_cov_event_date <= max_index_date & earliest_cov_event_date >= min_index_date ~ 1,
#     TRUE ~ 0
#   )) %>% 
#   mutate(had_covid_in_window = as.character(had_covid_in_window)) %>% 
#   ggplot() +
#   geom_density(aes(x=rsv_evidence_date, fill = had_covid_in_window), position="stack")

all_rsv_with_flu <- results_tbl("all_rsv_cases") %>% 
  group_by(person_id) %>% 
  slice_min(rsv_evidence_date, with_ties = FALSE) %>% 
  ungroup() %>% 
  mutate(min_index_date = as.Date(rsv_evidence_date - days(180))) %>% 
  mutate(max_index_date = as.Date(rsv_evidence_date - days(30))) %>% 
  get_flu_evidence_NN(min_date = cohort_entry_start_date, max_date = cohort_entry_end_date) %>% 
  compute_new(indexes=c("person_id"))
all_rsv_with_flu %>% output_tbl("all_rsv_w_flu")

full_disease_panel <- results_tbl("all_rsv_w_covid") %>% 
  select(person_id, rsv_evidence_date, min_index_date, max_index_date, earliest_cov_event_date, covid_event_positive) %>% 
  left_join(results_tbl("all_rsv_w_flu"), by="person_id") %>% 
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




# full_disease_panel %>% 
#   ggplot() +
#   geom_density(aes(x=rsv_evidence_date, fill=pre_acute_disease), position="stack")

# table1(~ pre_acute_disease, data = full_disease_panel)

## Make figure here for the table of the full landscape of infections
# Use the existing cohort
# Use the new RSV disease panel
# Join
# Make a figure

### NHAT HERE IS THE CODE FOR THE PRETTY FIGURE THAT SUCHITRA WANTS!

data <- results_tbl("resp_study_cohort_analytic_dataset") %>% 
  select(person_id, infection_date = ce_date, disease = sub_cohort) %>% 
  mutate(disease = if_else(disease == "COVID", "SARS-CoV-2", disease)) %>%
  union(results_tbl("descriptive_rsv_surge_data") %>% 
          select(person_id, infection_date = rsv_evidence_date) %>% 
          mutate(disease = "RSV")) %>%  collect() %>%
  mutate(disease = factor(disease, levels = c("SARS-CoV-2", "Influenza", "RSV")))%>%
  mutate(infection_week = floor_date(infection_date, unit="week")) %>%
  mutate(infection_week_ct = as.factor(difftime(infection_week, min(infection_week, na.rm = TRUE), unit = "weeks")))%>%
  group_by(infection_week, infection_week_ct, disease) %>% 
  summarise(n=n_distinct(person_id)) %>% 
  ungroup() %>% arrange(infection_week)

xtick_breaks <- c(as.numeric(difftime(as.Date("2021-07-01"),as.Date("2021-05-31"),unit = "weeks")),
                  as.numeric(difftime(as.Date("2021-10-01"),as.Date("2021-05-31"),unit = "weeks")),
                  as.numeric(difftime(as.Date("2022-01-01"),as.Date("2021-05-31"),unit = "weeks")),
                  as.numeric(difftime(as.Date("2022-04-01"),as.Date("2021-05-31"),unit = "weeks")),
                  as.numeric(difftime(as.Date("2022-07-01"),as.Date("2021-05-31"),unit = "weeks")),
                  as.numeric(difftime(as.Date("2022-10-01"),as.Date("2021-05-31"),unit = "weeks")),
                  as.numeric(difftime(as.Date("2023-01-01"),as.Date("2021-05-31"),unit = "weeks")),
                  as.numeric(difftime(as.Date("2023-04-01"),as.Date("2021-05-31"),unit = "weeks")))

xtick_labels <- c("2021-07", "2021-10", 
                  "2022-01", "2022-04", "2022-07", "2022-10",
                  "2023-01", "2023-04")

data %>% 
  ggplot() +
  geom_bar(aes(x=infection_week_ct, y=n, fill=disease), position="stack", stat="identity")+
  geom_vline(xintercept = as.numeric(difftime(as.Date("2022-03-01"),as.Date("2021-05-31"),unit = "weeks")), linetype = "dashed", color = "red")+
  geom_vline(xintercept = as.numeric(difftime(as.Date("2022-07-01"),as.Date("2021-05-31"),unit = "weeks")), linetype = "dashed", color = "red")+
  geom_vline(xintercept = factor(round(as.numeric(difftime(as.Date("2021-06-01"),as.Date("2021-05-31"),unit = "weeks")))), linetype = "dashed", color = "blue")+
  geom_vline(xintercept = as.numeric(difftime(as.Date("2022-10-01"),as.Date("2021-05-31"),unit = "weeks")), linetype = "dashed", color = "blue")+
  scale_fill_brewer(palette = "Dark2") +
  scale_x_discrete(breaks = factor(round(xtick_breaks)), labels = xtick_labels, limits = factor(-5:90)) +
  labs(x = "Infection week", y = "Counts (n)") + 
  theme_bw()

ggsave("figures/Figure_1_epi_curve.png", width = 10, height=6.5, "units" = "in")  

# analytic_dataset_resp %>% filter(sex_cat == "Other/unknown"| (!is.na(covid_index_date_imputed)& covid_index_date_imputed!=0)) %>%
# distinct_ct()





