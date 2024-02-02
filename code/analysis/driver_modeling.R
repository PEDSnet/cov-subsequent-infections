### Driver file for taking up an analytic dataset with outcomes and covariates, performing IPW weighting and Cox modeling
## for results

library(survey)
library(survival)
library(WeightIt)

cohort_outcomes <- results_tbl("test_cohort_outcomes") %>% 
  filter(sub_cohort %in% c("Covid", "Influenza")) %>% 
  mutate(exposure = case_when(sub_cohort=="Covid" ~ 1,
                              sub_cohort == "Influenza" ~ 0)) %>% 
  mutate(ce_month = as.Date(floor_date(ce_date, unit="month"))) %>% 
  mutate(ce_days_secular = as.numeric(ce_month - as.Date("2021-06-01"))) %>% 
  mutate(days_til_earliest_resp_outcome = ifelse(is.na(days_til_earliest_resp_outcome), 180, days_til_earliest_resp_outcome)) %>% 
  mutate(days_til_earliest_general_outcome = ifelse(is.na(days_til_earliest_general_outcome), 180, days_til_earliest_general_outcome)) %>% 
  mutate(days_til_earliest_rsv_outcome = ifelse(is.na(days_til_earliest_rsv_outcome), 180, days_til_earliest_rsv_outcome)) %>% 
  compute_new() %>% 
  collect()

weights.out <- weightit(as.formula("exposure~ce_days_secular"), data=cohort_outcomes,
                        method="glm", estimand="ATE") #focal = "circle pill")
cohort_outcomes$w = weights.out$weights

cohort_outcomes %>% 
  group_by(sub_cohort, has_postacute_general_outcome) %>% 
  summarise(n=n_distinct(person_id), weighted_n = sum(w)) %>% 
  ungroup() %>% 
  group_by(sub_cohort) %>% 
  mutate(cohort_total_w=sum(weighted_n)) %>% 
  mutate(prop_w=weighted_n/cohort_total_w) %>% 
  mutate(cohort_total=sum(n)) %>% 
  mutate(prop=n/cohort_total)

## TODO have an IPW plot here

cox_mod_adj <- coxph(Surv(days_til_earliest_general_outcome, has_postacute_general_outcome) ~ exposure + ce_days_secular, #need the time to event TODO start here on Tuesday
                     data=cohort_outcomes,
                     weights = w)

summary(cox_mod_adj)


