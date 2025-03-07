###### Anchor dates for conditions
postacute_min_start_date = as.Date("2022-02-01")
postacute_max_end_date = as.Date("2023-01-01")

cohort_entry_start_date = as.Date("2022-01-01")
cohort_entry_end_date = as.Date("2022-07-01")


### Step 1: Analysis of time-to-outcome for RSV outcomes

## Saved table with exclusion reasons

## TODO: analysis of test panel results by site, in cohort entry period, how many positive/negative pairs of tests are there (i.e. covid, influenza, RSV same day)
description = "covid_flu_nonexclude"

cohort_group <- results_tbl( paste0(description,'_tabc')) 

# flu_tests <- 
#   load_codeset("lab_influenza_filtered") %>% 
#   left_join(cdm_tbl("measurement_labs") %>% 
#               filter(measurement_date >= cohort_entry_start_date, measurement_date < postacute_max_end_date) %>% 
#               inner_join(cohort_group %>% select(person_id), by="person_id") %>% 
#               select(measurement_concept_id, measurement_source_concept_id, person_id, site, measurement_date, value_as_concept_id, value_source_value, value_as_number),
#             by=c("concept_id"="measurement_source_concept_id")) %>% 
#   compute_new(indexes=c("person_id"))
# 
# rsv_tests <- 
#   load_codeset("lab_rsv_final") %>% 
#   left_join(cdm_tbl("measurement_labs") %>% 
#               filter(measurement_date >= cohort_entry_start_date, measurement_date < postacute_max_end_date) %>% 
#               inner_join(cohort_group %>% select(person_id), by="person_id") %>% 
#               select(measurement_concept_id, person_id, site, measurement_date, value_as_concept_id),
#             by=c("concept_id"="measurement_concept_id")) %>% 
#   compute_new(indexes=c("person_id"))
# 
# covid_tests <- 
#   cdm_tbl("observation_derivation_recover") %>% 
#   inner_join(cohort_group %>% select(person_id), by="person_id") %>% 
#   filter(observation_date >= cohort_entry_start_date, observation_date < postacute_max_end_date) %>% 
#   filter(observation_concept_id %in% c(2000001530, 2000001529)) %>% 
#   select(person_id, site, measurement_date = observation_date, measurement_concept_id=observation_concept_id, 
#          value_as_concept_id) %>% 
#   group_by(person_id) %>%
#   compute_new(indexes=c("person_id"))
# 
# ####
# 
# ## NExt step: can combine peoples tests and then filter to just the day 0 patients to figure out what's going on here
# 
# ### Combined person, test, date, value
# combined_tests <-
#   flu_tests %>% 
#   group_by(person_id, site, concept_id, concept_name, concept_code, measurement_date, value_as_concept_id) %>% 
#   filter(row_number() == 1) %>% 
#   ungroup() %>% 
#   select(person_id, site, measurement_date, flu_concept_id = concept_id, flu_concept_name = concept_name, flu_concept_code = concept_code,
#          flu_value_as_concept_id = value_as_concept_id) %>% 
#   full_join(
#     covid_tests %>% 
#       group_by(person_id, site, measurement_concept_id, measurement_date, value_as_concept_id) %>% 
#       filter(row_number() == 1) %>% 
#       ungroup() %>% 
#       select(person_id, site, measurement_date, covid_concept_id = measurement_concept_id,
#              covid_value_as_concept_id = value_as_concept_id),
#     by=c("person_id", "site", "measurement_date")
#   ) %>% 
#   full_join(
#     rsv_tests %>% 
#       group_by(person_id, site, concept_id, concept_name, concept_code, measurement_date, value_as_concept_id) %>% 
#       filter(row_number() == 1) %>% 
#       ungroup() %>% 
#       select(person_id, site, measurement_date, rsv_concept_id = concept_id, rsv_concept_name = concept_name, rsv_concept_code = concept_code,
#              rsv_value_as_concept_id = value_as_concept_id),
#     by=c("person_id", "site", "measurement_date")
#   ) %>% 
#   compute_new(indexes=c("person_id", "site"))
# 
# combined_tests_collect <-
#   combined_tests %>% 
#   collect() %>% 
#   select(person_id, site, measurement_date,
#          flu_value_as_concept_id, covid_value_as_concept_id, rsv_value_as_concept_id) %>% 
#   # pivot_longer(cols = c("flu_concept_id", "covid_concept_id", "rsv_concept_id"), names_to = "virus_concept", values_to = "concept_id") %>% 
#   pivot_longer(cols = c("flu_value_as_concept_id", "covid_value_as_concept_id", "rsv_value_as_concept_id"),
#                names_to = "virus", values_to = "value_as_concept_id")
# 
# combined_tests_copy <- copy_to_new(df = combined_tests_collect)
# 
# combined_tests_labeled <-
#   combined_tests_copy %>% 
#   left_join(vocabulary_tbl("concept") %>% 
#               select(concept_id, value_as_concept_name = concept_name), by=c("value_as_concept_id"="concept_id")) %>% 
#   compute_new(indexes=c("person_id", "site"))
# 
# combined_tests_labeled %>% 
#   filter(!is.na(value_as_concept_id)) %>% 
#   group_by(site, person_id, measurement_date) %>% 
#   summarise(n_different_tests_on_date = n_distinct(virus)) %>% 
#   ungroup() %>% 
#   group_by(n_different_tests_on_date, site) %>% 
#   summarise(n_pats = n_distinct(person_id)) %>% 
#   ungroup() %>% 
#   ggplot() +
#   geom_bar(aes(x=n_different_tests_on_date, y = n_pats, fill=site), stat="identity", position="dodge") +
#   facet_wrap(~site)
# 
# ## Thing that we want to know: how many patients have a combo of poitive/negative for another thing
# combined_tests_labeled %>% 
#   filter(!is.na(value_as_concept_id)) %>% 
#   group_by(site, virus, value_as_concept_name, measurement_date) %>% 
#   summarise(n_pats = n_distinct(person_id))
# 
# combined_tests_labeled %>% 
#   filter(!is.na(value_as_concept_id)) %>% 
#   select(site, person_id, measurement_date, virus, value_as_concept_name) %>% 
#   pivot_wider(names_from = virus, values_from = value_as_concept_name) %>% 
#   group_by(site, measurement_date, covid_value_as_concept_id, flu_value_as_concept_id, rsv_value_as_concept_id) %>% 
#   summarise(n_pats = n_distinct(person_id))
# 
# 
# 
# test_panels <-
#   combined_tests_labeled %>% 
#   filter(!is.na(value_as_concept_id)) %>% 
#   mutate(test_result = case_when(
#     value_as_concept_name %in% c("No matching concept", "Unable", "Other",
#                                  "+", "33.3") ~ "NA or other",
#     value_as_concept_name %in% c("Detected", "Positive", "Equivocal", "Normal range",
#                                  "First positive instance") ~ "Positive",
#     value_as_concept_name %in% c("Not detected", "Negative", "COVID-19 Test Result Unknown", "Indeterminate", "Inconclusive") ~ "Negative/Unknown",
#     TRUE ~ "NA"
#   )) %>% 
#   select(site, person_id, measurement_date, virus, test_result) %>% 
#   group_by(site, person_id, measurement_date, virus, test_result) %>% 
#   filter(row_number()==1) %>% 
#   ungroup() %>% 
#   pivot_wider(id_cols = c(person_id, measurement_date, site),
#     names_from = virus, values_from = test_result, values_fill = "No test evidence") %>% 
#   compute_new() %>% 
#   mutate(panel_value = paste0("Cov: ", covid_value_as_concept_id, ", Flu: ", flu_value_as_concept_id, ", RSV: ", rsv_value_as_concept_id))
  
test_panels %>% output_tbl("cohort_test_panels", indexes=c("person_id"))

#### Would be great to have some type of model that classifies patients as what order infection they have from the data, and the probability

test_panels %>% 
  left_join(cohort_group %>% filter(lab_confirmed_rsv==1, days_from_index_to_rsv==0) %>%
              select(person_id, ce_date) %>% 
              mutate(pid=row_number()), by="person_id") %>% 
  mutate(days_between = ((as.numeric(measurement_date - ce_date)))) %>% 
  group_by(person_id) %>% 
  mutate(max_days_between = max(abs(days_between))) %>% 
  ungroup() %>% 
  mutate(max_group = max_days_between %% 4) %>% 
  filter(abs(days_between) <= 30) %>% ggplot() +
  geom_line(aes(x=days_between, y = pid, group=pid))  +
  geom_point(aes(x=days_between, y = pid, group=pid, color=panel_value))  +
  facet_wrap(~max_group, scales="free")

cohort_group %>% filter(lab_confirmed_rsv==1, days_from_index_to_rsv==0) %>%
  group_by(is.na(cov_test_date), is.na(covid_dx_date), covid_index_date_imputed) %>% 
  summarise(n=n_distinct(person_id)) %>% 
  ungroup() %>% 
  mutate(total=sum(n)) %>% 
  mutate(n/total)

## How to triage patients real covid date as compared to rsv



####

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
  output_tbl("test_panel_aws_summaries_cohort")

p <- panel_summaries %>% 
  filter(!panel_value=="Cov: No test evidence, Flu: No test evidence, RSV: No test evidence") %>% 
  # filter(site=="chop") %>% 
  # filter(person_id == 53606) %>% 
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
  facet_wrap(~site, ncol = 5) +
  theme_bw() 

leg <- get_legend(p)  
p + theme(legend.position = "None")
as_ggplot(leg)
ggsave("results/test_panel_summary_cohort.jpg", width = 20, height = 12, dpi=800)


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
  group_by(site, covid_value_as_concept_name, flu_value_as_concept_name, rsv_value_as_concept_name)) 
  summarise(n_distinct(person_id)) %>% 
  ggplot() +
  geom_bar(aes(x=))



