library(survey)
library(survival)
library(WeightIt)
library(cobalt)
library(ggplot2)


postacute_min_start_date = as.Date("2022-02-01")
postacute_max_end_date = as.Date("2023-01-01")

cohort_entry_start_date = as.Date("2022-01-01")
cohort_entry_end_date = as.Date("2022-07-01")

cohort_1_label = "_rsv_study_cohort_"
cohort_1_label = "sample_pats"

exclusion_reasons <- results_tbl(paste0(cohort_1_label, "_w_exclusion_reasons"))

attrition <- results_tbl(paste0(cohort_1_label, "_attrition"))

final_cohort_demo <- results_tbl(paste0(cohort_1_label, "_cohort_demo"))

rsv_outcomes <- results_tbl(paste0(cohort_1_label, "_rsv_outcomes"))

cohort_with_rsv_outcomes <-
  final_cohort_demo %>% 
  left_join(rsv_outcomes, by="person_id")

analytic_dataset <-
  results_tbl(paste0(cohort_1_label, "_analytic_dataset"))


## TODO why are there multiple rows for the same person
# TODO the RSV function is saving rows for each RSV index date. That makes sense, because need to track all infections, so might need some fancy logic for exlusions though
# not just earliest because there could be like a prior case, and then a pot-acute case
final_analytic_dataset <-
  analytic_dataset %>% 
  filter(is.na(covid_index_date_imputed) | covid_index_date_imputed==0) %>% 
  group_by(person_id) %>% 
  mutate(exclude = max(exclude_based_on_rsv_outcome)) %>% 
  filter(exclude < 1) %>% 
  slice_min(rsv_evidence_date, with_ties=FALSE) %>% 
  ungroup() %>% 
  compute_new(indexes=c("person_id"))

#### IPTW

iptw_dataset <-
  final_analytic_dataset %>% 
  mutate(ce_week = as.Date(floor_date(ce_date, unit="week"))) %>% 
  mutate(study_start_date = as.Date(cohort_entry_start_date)) %>% 
  mutate(ce_date_days_numeric = as.numeric(ce_date-study_start_date)) %>% 
  mutate(ce_week_days_numeric = as.numeric(ce_week-study_start_date)) %>% 
  compute_new(indexes=c("person_id"))

### Now first want to limit the cohort to the COVID/flu comparison:
iptw_covid_flu <-
  iptw_dataset %>% 
  filter(sub_cohort %in% c("COVID", "Influenza")) %>% 
  mutate(covid = ifelse(sub_cohort=="COVID", 1, 0)) %>% 
  mutate(rsv_outcome = case_when(
    is.na(rsv_in_post_acute_period) ~ 0,
    rsv_in_post_acute_period == "No" ~ 0,
    rsv_in_post_acute_period == "Yes" ~ 1
  )) %>% 
  collect()


#### Calling modeling stuff from function
weight_formula = as.formula('covid~age_group+sex_cat+race_eth_cat+ce_week_days_numeric')

covid_flu_iptw_results <- 
  make_iptw_and_mw_weights_colby(weight_formula = weight_formula,
                               dataset = iptw_covid_flu,
                               description = "covid_flu_sample"
                               )


#iptw_covid_flu %>% group_by(ce_week_days_numeric, sub_cohort) %>% summarise(n=n()) %>% ggplot() + geom_bar(aes(x=ce_week_days_numeric, y = n, fill=sub_cohort), stat="identity", position="dodge")

# bal.tab(exposure ~ ce_week_days_numeric + sex_cat + race_eth_cat,
#         data = iptw_covid_flu, estimand = "ATT", thresholds = c(m = .05))
# 
# weights_result <- weightit(as.formula("exposure~ce_week_days_numeric + sex_cat + race_eth_cat"), 
#                            data=iptw_covid_flu,
#                         method="glm", estimand="ATT") #focal = "circle pill")
# 
# iptw_covid_flu$wt = weights_result$weights
# 
# summary(weights_result)
# 
# bal.tab(weights_result, stats = c("m", "v"), thresholds = c(m = .05))
# 
# weighted_entropy <- weightit(exposure ~ ce_week_days_numeric + sex_cat + race_eth_cat,
#                   data = iptw_covid_flu, estimand = "ATT", method = "ebal")
# summary(weighted_entropy)
# 
# bal.tab(weighted_entropy, stats = c("m", "v"), thresholds = c(m = .05))
# 
# ## Plain IPTW proportion outcome:
# 
# iptw_covid_flu %>% 
#   group_by(sub_cohort, rsv_outcome) %>% 
#   summarise(n=n_distinct(person_id), weighted_n = sum(wt)) %>% 
#   ungroup() %>% 
#   group_by(sub_cohort) %>% 
#   mutate(cohort_total_w=sum(weighted_n)) %>% 
#   mutate(prop_w=weighted_n/cohort_total_w) %>% 
#   mutate(cohort_total=sum(n)) %>% 
#   mutate(prop=n/cohort_total)

## Checking balance




#### Logit modeling
## Create a model for the rsv outcome using the weighted data and the covariates
## Create an unweighted logistic regression for rsv outcome using unweighted data and the covariates
## Create an unweighted model, unadjusted with no covariates
logit_no_adjustment <- glm(rsv_outcome ~ covid, 
               data = covid_flu_iptw_results$dataset, family = quasibinomial(link="logit"), 
               weights = capped_iptw)

summary(logit_no_adjustment)

logit_adjustment <- glm(rsv_outcome ~ covid + race_eth_cat + age_group + sex_cat, 
                           data = covid_flu_iptw_results$dataset, family = quasibinomial(link="logit"), 
                           weights = capped_iptw)

summary(logit_adjustment)


#### Cox regression






