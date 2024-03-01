###### Anchor dates for conditions
postacute_min_start_date = as.Date("2022-04-01")
postacute_max_end_date = as.Date("2023-07-01")

cohort_entry_start_date = as.Date("2022-03-01")
cohort_entry_end_date = as.Date("2023-01-01")


### Step 1: Analysis of time-to-outcome for RSV outcomes

## Saved table with exclusion reasons

## TODO: analysis of test panel results by site, in cohort entry period, how many positive/negative pairs of tests are there (i.e. covid, influenza, RSV same day)

flu_tests <- 
  load_codeset("lab_influenza_filtered") %>% 
  left_join(cdm_tbl("measurement_labs") %>% 
              filter(measurement_date >= cohort_entry_start_date, measurement_date < postacute_max_end_date) %>% 
              select(measurement_concept_id, person_id, site, measurement_date, value_as_concept_id),
            by=c("concept_id"="measurement_concept_id")) %>% 
  compute_new(indexes=c("person_id"))

rsv_tests <- 
  load_codeset("lab_rsv_final") %>% 
  left_join(cdm_tbl("measurement_labs") %>% 
              filter(measurement_date >= cohort_entry_start_date, measurement_date < postacute_max_end_date) %>% 
              select(measurement_concept_id, person_id, site, measurement_date, value_as_concept_id),
            by=c("concept_id"="measurement_concept_id")) %>% 
  compute_new(indexes=c("person_id"))

covid_tests <- 
  cdm_tbl("observation_derivation_recover") %>% 
  filter(observation_date >= cohort_entry_start_date, observation_date < postacute_max_end_date) %>% 
  filter(observation_concept_id %in% c(2000001530, 2000001529)) %>% 
  select(person_id, site, measurement_date = observation_date, measurement_concept_id=observation_concept_id, 
         value_as_concept_id) %>% 
  group_by(person_id) %>%
  compute_new(indexes=c("person_id"))

####

### Combined person, test, date, value
combined_tests <-
  flu_tests %>% 
  distinct(person_id, site, concept_id, concept_name, concept_code, measurement_date, value_as_concept_id) %>% 
  select(person_id, site, measurement_date, flu_concept_id = concept_id, flu_concept_name = concept_name, flu_concept_code = concept_code,
         flu_value_as_concept_id = value_as_concept_id) %>% 
  full_join(
    covid_tests %>% 
      distinct(person_id, site, measurement_concept_id, measurement_date, value_as_concept_id) %>% 
      select(person_id, site, measurement_date, covid_concept_id = measurement_concept_id,
             covid_value_as_concept_id = value_as_concept_id),
    by=c("person_id", "site", "measurement_date")
  ) %>% 
  full_join(
    rsv_tests %>% 
      distinct(person_id, site, concept_id, concept_name, concept_code, measurement_date, value_as_concept_id) %>% 
      select(person_id, site, measurement_date, rsv_concept_id = concept_id, rsv_concept_name = concept_name, rsv_concept_code = concept_code,
             rsv_value_as_concept_id = value_as_concept_id),
    by=c("person_id", "site", "measurement_date")
  ) %>% 
  compute_new(indexes=c("person_id", "site"))

combined_tests_labeled <-
  combined_tests %>% 
  select(person_id, site, measurement_date,
         flu_value_as_concept_id, covid_value_as_concept_id, rsv_value_as_concept_id) %>% 
  # pivot_longer(cols = c("flu_concept_id", "covid_concept_id", "rsv_concept_id"), names_to = "virus_concept", values_to = "concept_id") %>% 
  pivot_longer(cols = c("flu_value_as_concept_id", "covid_value_as_concept_id", "rsv_value_as_concept_id"), names_to = "virus", values_to = "value_as_concept_id") %>% 
  left_join(vocabulary_tbl("concept") %>% 
              select(concept_id, value_as_concept_name = concept_name), by=c("value_as_concept_id"="concept_id")) %>% 
  compute_new(indexes=c("person_id", "site"))

combined_tests_labeled %>% 
  filter(!is.na(value_as_concept_id)) %>% 
  group_by(site, person_id, measurement_date) %>% 
  summarise(n_different_tests_on_date = n_distinct(virus)) %>% 
  ungroup() %>% 
  group_by(n_different_tests_on_date, site) %>% 
  summarise(n_pats = n_distinct(person_id)) %>% 
  ungroup() %>% 
  ggplot() +
  geom_bar(aes(x=n_different_tests_on_date, y = n_pats, fill=site), stat="identity", position="dodge") +
  facet_wrap(~site)

## Thing that we want to know: how many patients have a combo of poitive/negative for another thing
combined_tests_labeled %>% 
  filter(!is.na(value_as_concept_id)) %>% 
  group_by(site, virus, value_as_concept_name, measurement_date) %>% 
  summarise(n_pats = n_distinct(person_id))

combined_tests_labeled %>% 
  filter(!is.na(value_as_concept_id)) %>% 
  select(site, person_id, measurement_date, virus, value_as_concept_name) %>% 
  pivot_wider(names_from = virus, values_from = value_as_concept_name) %>% 
  group_by(site, measurement_date, covid_value_as_concept_id, flu_value_as_concept_id, rsv_value_as_concept_id) %>% 
  summarise(n_pats = n_distinct(person_id))



test_panels <-
  combined_tests_labeled %>% 
  filter(!is.na(value_as_concept_id)) %>% 
  mutate(test_result = case_when(
    value_as_concept_name %in% c("No matching concept", "Unable", "Other",
                                 "+", "33.3") ~ "NA or other",
    value_as_concept_name %in% c("Detected", "Positive", "Equivocal", "Normal range",
                                 "First positive instance") ~ "Positive",
    value_as_concept_name %in% c("Not detected", "Negative", "COVID-19 Test Result Unknown", "Indeterminate", "Inconclusive") ~ "Negative/Unknown",
    TRUE ~ "NA"
  )) %>% 
  select(site, person_id, measurement_date, virus, test_result) %>% 
  pivot_wider(names_from = virus, values_from = test_result, values_fill = "No test evidence") %>% 
  compute_new() %>% 
  mutate(panel_value = paste0("Cov: ", covid_value_as_concept_id, ", Flu: ", flu_value_as_concept_id, ", RSV: ", rsv_value_as_concept_id))
  
panel_summaries <-
  test_panels %>% 
  group_by(site, panel_value, measurement_date) %>% 
  summarise(n_pats = n_distinct(person_id)) %>% 
  ungroup() %>% 
  group_by(site, measurement_date) %>% 
  mutate(n_panels_on_date = sum(n_pats)) %>% 
  ungroup() %>% 
  mutate(fraction_with_value = n_pats/n_panels_on_date) %>% 
  compute_new()

panel_summaries %>% 
  output_tbl("test_panel_pedsnet_summaries")

panel_summaries %>% 
  filter(!panel_value=="Cov: No test evidence, Flu: No test evidence, RSV: No test evidence") %>% 
  # filter(panel_value %in% c(
  #   "Cov: Positive, Flu: Positive, RSV: Positive",
  #   "Cov: Positive, Flu: Positive, RSV: No test evidence",
  #   "Cov: No test evidence, Flu: Positive, RSV: No test evidence",
  #   "Cov: No test evidence, Flu: Positive, RSV: Positive",
  #   "Cov: No test evidence, Flu: No test evidence, RSV: Positive",
  #   "Cov: Positive, Flu: No test evidence, RSV: Positive",
  #   "Cov: Negative, Flu: No test evidence, RSV: No test evidence",
  #   "Cov: No test evidence, Flu: Negative, RSV: No test evidence",
  #   "Cov: No test evidence, Flu: No test evidence, RSV: Negative"
  # )) %>% 
  ggplot() +
  geom_bar(aes(x=measurement_date, y = fraction_with_value, fill = panel_value), 
           stat="identity", 
           position="stack") +
  facet_wrap(~site) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("results/test_panel_summary.jpg", width = 20, height = 12, dpi=500)


### think about actual plot we want:
# Per site, what are the proportions of patients with each panel result (same day) out of all the tests that day
# And then what are the proportions over all of panel results 




# combined_tests_labeled <-
#   combined_tests %>% 
#   left_join(vocabulary_tbl("concept") %>% 
#               select(concept_id, covid_value_as_concept_name = concept_name),
#             by=c("covid_value_as_concept_id"="concept_id")) %>% 
#   left_join(vocabulary_tbl("concept") %>% 
#               select(concept_id, flu_value_as_concept_name = concept_name),
#             by=c("flu_value_as_concept_id"="concept_id")) %>% 
#   left_join(vocabulary_tbl("concept") %>% 
#               select(concept_id, rsv_value_as_concept_name = concept_name),
#             by=c("rsv_value_as_concept_id"="concept_id")) %>%
#   compute_new(indexes=c("person_id", "site"))
#   

### Plots to make:
## For tests on the same day, what proportion are positive/negative for the different things

## First: for each patient, how many tests a day do they have:
combined_tests_labeled %>% 
  group_by(person_id) %>% 
  summarise(n_)

combined_tests_labeled %>% 
  group_by(site, covid_value_as_concept_name, flu_value_as_concept_name, rsv_value_as_concept_name) %>% 
  summarise(n_distinct(person_id)) %>% 
  ggplot() +
  geom_bar(aes(x=))


