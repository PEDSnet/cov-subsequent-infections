### File for constructing code to generate variables given a cohort
## Includes generation of intermediate tables but not functions for permanent output, which should be computed separately.


## Variable functions that need to be made:
variables <- c("rsv_diagnosis", "rsv_lab_positive", 
               "respiratory_lab_positive", "respiratory_infection_diagnosis",
               "any_infection_lab_positive", "any_infection_diagnosis",
               "covid_diagnosis", "covid_lab_positive",
               "influenza_diagnosis", "influenza_lab_positive",
               "age", "sex", "race_ethnicity", "medical_complexity", 
               "severity", "insurance", "insurance_carrier",
               "utilization", "vaccination")


## Analytic dataset

# Format design: Tidy data set. 
## Dataset ultimately needs to be formatted for a Cox model, but also a logistic regression and a Poisson regression, though those can
# use functions like group by and summarise to group by person and summarise counts, so that can actually come later. The Cox df will be the main one

## Cox DF spec:
## 
cox_df_cols <- c("person_id", "ce_date", "sub_cohort", #base cohort info
                 "covariate", # co-variates
                 "respiratory_outcome", "general_infection_outcome", "rsv_outcome", # outcomes
                 "days_to_respiratory_outcome", "days_to_general_infection_outcome", "days_to_rsv_outcome", # time to outcomes
                 )


## Write code to apply each of these variables to the cohort

## Getting outcomes:
## respiratory infection
base_cohort <- results_tbl("base_cohort") %>% 
  group_by(sub_cohort) %>% 
  slice_sample(n=100) %>% 
  ungroup() %>% 
  compute_new()

base_cohort_resp_outcomes_test <-
  base_cohort %>% 
  flag_outcome_resp_infections(outcome_start_date = as.Date("2021-07-01"),
                               outcome_end_date = as.Date("2023-08-01"))

base_cohort_any_outcomes_test <-
  base_cohort %>% 
  flag_outcome_infections(outcome_start_date = as.Date("2021-07-01"),
                               outcome_end_date = as.Date("2023-08-01"))

base_cohort_rsv_outcomes_test <-
  base_cohort %>% 
  flag_rsv_outcome_draft(outcome_start_date = as.Date("2022-04-01"),
                               outcome_end_date = as.Date("2023-07-01"))

### Then need to edit for records for patients within the proper cohort entry period

## First, EDA. How many infections does each person have?
base_cohort_resp_outcomes_test %>% 
  group_by(person_id, condition_start_date) %>% 
  filter(row_number()==1) %>% 
  ungroup() %>% 
  group_by(person_id) %>% 
  summarise(unique_dates_of_infection = n()) %>% 
  ggplot() +
  geom_histogram(aes(x=unique_dates_of_infection))

base_cohort_rsv_outcomes_test %>% 
  group_by(person_id, earliest_rsv_evidence) %>% 
  filter(row_number()==1) %>% 
  ungroup() %>% 
  group_by(person_id) %>% 
  summarise(unique_dates_of_infection = n()) %>% 
  ggplot() +
  geom_histogram(aes(x=unique_dates_of_infection))


base_cohort_any_outcomes_test %>% 
  group_by(person_id, condition_start_date) %>% 
  filter(row_number()==1) %>% 
  ungroup() %>% 
  group_by(person_id) %>% 
  summarise(unique_dates_of_infection = n()) %>% 
  ggplot() +
  geom_histogram(aes(x=unique_dates_of_infection))

### Next, take first evidence of outcome, since most people have just 1, and then half have 2 or more.
## TODO note that the earliest infection here is due to influenza on the same day, for that cohort. 
## So, need to make it such that we consider the earliest infection IN the post-acute period?
## Would be helpful to make that trajectory plot of all diagnosis evidence, though, because there could be follow-ups & we need a gap between index infection and next infection
earliest_resp_infection <-
  base_cohort_resp_outcomes_test %>% 
  left_join(base_cohort, by="person_id") %>% 
  mutate(days_til_resp_outcome = as.numeric(condition_start_date - ce_date)) %>% 
  filter(days_til_resp_outcome > 7) %>% 
  group_by(person_id) %>% 
  slice_min(condition_start_date, with_ties = FALSE) %>% 
  select(person_id, earliest_resp_outcome = condition_start_date, resp_concept = concept_name) %>% 
  ungroup() %>% 
  compute_new()

cohort_with_resp_outcomes <-
  base_cohort %>% 
  left_join(earliest_resp_infection, by="person_id") %>% 
  mutate(days_til_earliest_resp_outcome = as.numeric(earliest_resp_outcome - ce_date)) %>% 
  compute_new()

cohort_with_resp_outcomes %>% 
  filter(days_til_earliest_resp_outcome >= 30, days_til_earliest_resp_outcome < 180) %>%
  ggplot() +
  geom_density(aes(x=days_til_earliest_resp_outcome, fill = sub_cohort), alpha=0.5)

cohort_with_resp_outcomes %>% 
  mutate(has_postacute_resp_outcome = case_when(
    !is.na(resp_concept) & days_til_earliest_resp_outcome >= 30 & days_til_earliest_resp_outcome < 180 ~ "Yes",
    TRUE ~ "No"
  )) %>% 
  filter(ce_date >= as.Date("2021-06-01"), ce_date < as.Date("2023-02-01")) %>% 
  group_by(has_postacute_resp_outcome, sub_cohort) %>% 
  summarise(n = n_distinct(person_id)) %>% 
  ggplot() +
  geom_bar(aes(x=sub_cohort, y = n, fill=has_postacute_resp_outcome), stat="identity", position="dodge")


earliest_general_infection <-
  base_cohort_any_outcomes_test %>% 
  left_join(base_cohort, by="person_id") %>% 
  mutate(days_til_any_outcome = as.numeric(condition_start_date - ce_date)) %>% 
  filter(days_til_any_outcome > 7) %>% 
  group_by(person_id) %>% 
  slice_min(condition_start_date, with_ties = FALSE) %>% 
  select(person_id, earliest_general_outcome = condition_start_date, general_concept = concept_name) %>% 
  ungroup() %>% 
  compute_new()

cohort_with_general_outcomes <-
  cohort_with_resp_outcomes %>% 
  left_join(earliest_general_infection, by="person_id") %>% 
  mutate(days_til_earliest_general_outcome = as.numeric(earliest_general_outcome - ce_date)) %>% 
  compute_new()

cohort_with_all_outcomes <-
  cohort_with_general_outcomes %>% 
  left_join(base_cohort_rsv_outcomes_test %>% 
              select(person_id, earliest_rsv_evidence), by = "person_id") %>% 
  mutate(days_til_earliest_rsv_outcome = as.numeric(earliest_rsv_evidence - ce_date)) %>% 
  mutate(has_postacute_resp_outcome = case_when(
    !is.na(resp_concept) & days_til_earliest_resp_outcome >= 30 & days_til_earliest_resp_outcome < 180 ~ 1,
    TRUE ~ 0
  )) %>% 
  mutate(has_postacute_general_outcome = case_when(
    !is.na(general_concept) & days_til_earliest_general_outcome >= 30 & days_til_earliest_general_outcome < 180 ~ 1,
    TRUE ~ 0
  )) %>% 
  mutate(has_postacute_rsv_outcome = case_when(
    !is.na(earliest_rsv_evidence) & days_til_earliest_rsv_outcome >= 30 & days_til_earliest_rsv_outcome < 180 ~ 1,
    TRUE ~ 0
  )) %>% 
  compute_new()

### time series of outcomes
cohort_with_all_outcomes %>% 
  filter(sub_cohort %in% c("Covid", "Influenza")) %>% 
  pivot_longer(cols = c(days_til_earliest_general_outcome, days_til_earliest_resp_outcome, days_til_earliest_rsv_outcome),
               names_to = "outcome", values_to = "days_until_outcome") %>% 
  ggplot() +
  geom_density(aes(x=days_until_outcome, fill=sub_cohort)) +
  facet_wrap(~outcome, scales = "free") +
  theme_bw()

## TODO this is a good intermiediary plot (could go in a "variable development" report, but shows we need some more logic around
## identifying outcome infections. Worried about the distribution being smushed at 30 days, it means, there's a lot behind that,
## that we are arbitrarily cutting off, which could be follow up or repeat diagnoses


## Next plot: time series by cohort entry date of proportion of people with the outcome, by cohort entry week
cohort_with_all_outcomes %>% 
  filter(sub_cohort %in% c("Covid", "Influenza")) %>% 
  pivot_longer(cols = c(days_til_earliest_general_outcome, days_til_earliest_resp_outcome, days_til_earliest_rsv_outcome),
               names_to = "outcome", values_to = "days_until_outcome") %>% 
  mutate(ce_month = floor_date(ce_date, unit = "month")) %>% 
  mutate(ce_week = floor_date(ce_date, unit = "week")) %>% 
  mutate(has_outcome = ifelse(is.na(days_until_outcome), 0, 1)) %>% 
  group_by(sub_cohort, ce_month, outcome, has_outcome) %>% 
  summarise(n = n_distinct(person_id)) %>% 
  ungroup() %>% 
  group_by(sub_cohort, ce_month, outcome) %>% 
  mutate(total_month = sum(n)) %>% 
  mutate(prop = as.numeric(n)/as.numeric(total_month)) %>% 
  filter(has_outcome==1) %>% 
  ggplot() +
  geom_point(aes(x=ce_month, y = prop, color = sub_cohort)) +
  geom_smooth(aes(x=ce_month, y = prop, color = sub_cohort), span = 0.3) +
  facet_wrap(~outcome) +
  theme_bw()


cohort_with_all_outcomes %>% 
  output_tbl("test_cohort_outcomes")
## USs this structure to build up cox model


## Then for the EDA, would be good to first look at DQ / availability based on each variable/codeset
## Then for analytical EDA, would be good to look at correlations and associations, pair-plots, etc. to see if there are associations between
# variables before putting them into a model, or finding interesting things to begin with. 

### Questions: Insurance, what visit do we want to use for `visit_payer` to figure out their insurance class?
## Seems like its this so far:
# plan_class        
# 1 Other public      
# 2 Private/Commercial
# 3 Other/Unknown     
# 4 Medicaid/sCHIP    
# 5 Medicare          
# 6 Self-pay  

## The beginnings of some insurance DQ:
results_tbl("base_cohort") %>% 
  left_join(cdm_tbl("visit_occurrence") %>% 
              select(visit_occurrence_id, person_id, visit_start_date), by=c("person_id", "ce_date"="visit_start_date")) %>% 
  left_join(cdm_tbl("visit_payer") %>% 
              select(visit_occurrence_id, plan_class), by="visit_occurrence_id") %>% 
  group_by(plan_class) %>% 
  summarise(n())




