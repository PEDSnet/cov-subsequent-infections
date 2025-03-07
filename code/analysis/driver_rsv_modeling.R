library(survey)
library(survival)
library(WeightIt)
library(cobalt)
library(ggplot2)


postacute_min_start_date = as.Date("2022-04-01")
postacute_max_end_date = as.Date("2023-01-01")

cohort_entry_start_date = as.Date("2022-03-01")
cohort_entry_end_date = as.Date("2022-07-01")

cohort_1_label = "rsv_study_cohort"
# cohort_1_label = "sample_pats"


analytic_dataset <-
  results_tbl(paste0(cohort_1_label, "_analytic_dataset")) %>% 
  filter(ce_date >= cohort_entry_start_date,
         ce_date < cohort_entry_end_date)


## TODO why are there multiple rows for the same person
# TODO the RSV function is saving rows for each RSV index date. That makes sense, because need to track all infections, so might need some fancy logic for exlusions though
# not just earliest because there could be like a prior case, and then a pot-acute case
final_analytic_dataset <-
  analytic_dataset %>% 
  filter(sex_cat != "Other/unknown") %>% 
  filter(is.na(covid_index_date_imputed) | covid_index_date_imputed==0) %>% 
  group_by(person_id) %>% 
  mutate(exclude = max(exclude_based_on_rsv_outcome)) %>% ## TODO just exclude for prior rsv
  filter(exclude < 1) %>% 
  slice_min(rsv_evidence_date, with_ties=FALSE) %>% 
  ungroup() %>% 
  # mutate(new_rsv_outcome_fu = case_when(
  #   rsv_occurrence_period %in% c("post acute (1-6 months after index)", "Future RSV 180-300 days)") ~ 1,
  #   TRUE ~ 0
  # )) %>% 
  compute_new(indexes=c("person_id"))

final_analytic_dataset_nonexclusion <-
  analytic_dataset %>% 
  filter(sex_cat != "Other/unknown") %>% 
  filter(is.na(covid_index_date_imputed) | covid_index_date_imputed==0) %>% 
  group_by(person_id) %>% 
  mutate(exclude = max(exclude_for_prior_rsv)) %>% 
  filter(exclude < 1) %>%
  slice_min(rsv_evidence_date, with_ties=FALSE) %>% 
  ungroup() %>% 
  # mutate(new_rsv_outcome_fu = case_when(
  #   rsv_occurrence_period %in% c("post acute (1-6 months after index)", "Future RSV 180-300 days)") ~ 1,
  #   TRUE ~ 0
  # )) %>% 
  compute_new(indexes=c("person_id"))

iptw_flu <- final_analytic_dataset_nonexclusion %>% 
  do_iptw(cohort_entry_start_date = cohort_entry_start_date,
          comparison_cohort = "Influenza",
          weight_formula = "covid~age_group+sex_cat+race_eth_cat+ce_week_factor",
          iptw_description = "covid_flu_comparison")

iptw_ari <- final_analytic_dataset_nonexclusion %>% 
  do_iptw(cohort_entry_start_date = cohort_entry_start_date,
          comparison_cohort = "Respiratory",
          weight_formula = "covid~age_group+sex_cat+race_eth_cat+ce_week_factor",
          iptw_description = "covid_ari_comparison")

#### IPTW

iptw_dataset <-
  final_analytic_dataset_nonexclusion %>% 
  mutate(ce_week = as.Date(floor_date(ce_date, unit="week"))) %>% 
  mutate(study_start_date = as.Date(cohort_entry_start_date)) %>% 
  mutate(ce_date_days_numeric = as.numeric(ce_date-study_start_date)) %>% 
  mutate(ce_week_days_numeric = as.numeric(ce_week-study_start_date)) %>% 
  collect() %>% 
  mutate(ce_week_factor = as.factor(ce_week_days_numeric)) 
  # compute_new(indexes=c("person_id"))

### Now first want to limit the cohort to the COVID/flu comparison:
iptw_covid_flu <-
  iptw_dataset %>% 
  filter(sub_cohort %in% c("COVID", "Influenza")) %>% 
  mutate(covid = ifelse(sub_cohort=="COVID", 1, 0)) 


#### Calling modeling stuff from function
weight_formula = as.formula('covid~age_group+sex_cat+race_eth_cat+ce_week_factor')

weight_formula_with_util = as.formula('covid~age_group+sex_cat+race_eth_cat+ce_week_factor+util_other+util_inpatient+util_outpatient+util_ed')

# TODO change the function also to output model coefficient results based on pre-weight analysis
covid_flu_iptw_results <- 
  make_iptw_and_mw_weights_colby(weight_formula = weight_formula_with_util,
                               dataset = iptw_covid_flu,
                               description = "covid_flu_nonexclude" # covid_flu_noutil
                               )
covid_flu_iptw_results$love_plot

ggsave(paste0("results/", "covid_flu_nonexclude", "_love_plot.png"), width = 8, height = 8, dpi = 400)

covid_flu_iptw_results <- 
  make_iptw_and_mw_weights_colby(weight_formula = weight_formula_with_util,
                                 dataset = iptw_covid_flu,
                                 description = "covid_flu_util"
  )


covid_flu_iptw_results$love_plot

ggsave(paste0("results/", "covid_flu_util", "_love_plot.png"), width = 8, height = 8, dpi = 400)


### Now first want to limit the cohort to the COVID/respiratory comparison:
iptw_covid_ari <-
  iptw_dataset %>% 
  filter(sub_cohort %in% c("COVID", "Respiratory")) %>% 
  mutate(covid = ifelse(sub_cohort=="COVID", 1, 0)) 


#### Calling modeling stuff from function
weight_formula = as.formula('covid~age_group+sex_cat+race_eth_cat+ce_week_factor')

weight_formula_with_util = as.formula('covid~age_group+sex_cat+race_eth_cat+ce_week_factor+util_other+util_inpatient+util_outpatient+util_ed')

# TODO change the function also to output model coefficient results based on pre-weight analysis
covid_ari_iptw_results <- 
  make_iptw_and_mw_weights_colby(weight_formula = weight_formula_with_util,
                                 dataset = iptw_covid_ari,
                                 description = "covid_ari_util"
  )


covid_ari_iptw_results$love_plot

ggsave(paste0("results/", "covid_ari_util", "_love_plot.png"), width = 8, height = 8, dpi = 400)



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

library(jtools)


models = list()

models$model_unweighted_unadjusted <- run_logistic_regression(data = covid_flu_iptw_results$dataset,
                                                       formula = as.formula("rsv_outcome ~ covid"),
                                                       use_weights = FALSE,
                                                       weight_col = NA)

models$model_weighted_lr_unadjusted <- run_logistic_regression(data = covid_flu_iptw_results$dataset,
                                                       formula = as.formula("rsv_outcome ~ covid"),
                                                       use_weights = TRUE,
                                                       weight_col = capped_iptw)

models$model_weighted_gbm_unadjusted <- run_logistic_regression(data = covid_flu_iptw_results$dataset,
                                                        formula = as.formula("rsv_outcome ~ covid"),
                                                        use_weights = TRUE,
                                                        weight_col = iptw_gbm)

models$model_weighted_lr_adjusted <- run_logistic_regression(data = covid_flu_iptw_results$dataset,
                                                         formula = as.formula("rsv_outcome ~ covid + race_eth_cat + age_group + sex_cat + ce_week_days_numeric"),
                                                         use_weights = TRUE,
                                                         weight_col = capped_iptw)

models$model_weighted_gbm_adjusted <- run_logistic_regression(data = covid_flu_iptw_results$dataset,
                                                         formula = as.formula("rsv_outcome ~ covid + race_eth_cat + age_group + sex_cat + ce_week_days_numeric"),
                                                         use_weights = TRUE,
                                                         weight_col = iptw_gbm)

plot_summs(models$model_unweighted_unadjusted, 
           models$model_weighted_lr_unadjusted,
           models$model_weighted_gbm_unadjusted,
           models$model_weighted_lr_adjusted,
           models$model_weighted_gbm_adjusted,
           model.names = c("Unweighted, unadjusted", "Weighted (LR), unadjusted", 
                           "Weighted (GBM), unadjusted", "Weighted (LR), adjusted",
                           "Weighted (GBM), adjusted"))

logit_adjustment <- glm(rsv_outcome ~ covid + race_eth_cat + age_group + sex_cat, 
                           data = covid_flu_iptw_results$dataset, family = quasibinomial(link="logit"), 
                           weights = capped_iptw)

summary(logit_adjustment)

logit_no_weighted <- glm(rsv_outcome ~ covid, 
                         data = covid_flu_iptw_results$dataset, family = binomial(link="logit"))

summary(logit_no_weighted)
## Unweighted, way more likely to get RSV after covid. But then weighted, with o adjustment, says less likely,
# but not statistically significant

logit_adjustment <- glm(rsv_outcome ~ covid + race_eth_cat + age_group + sex_cat + ce_week_days_numeric, 
                        data = covid_flu_iptw_results$dataset, family = quasibinomial(link="logit"), 
                        weights = capped_iptw)

summary(logit_adjustment)

logit_gbm <- glm(rsv_outcome ~ covid, 
                        data = covid_flu_iptw_results$dataset, family = quasibinomial(link="logit"), 
                        weights = iptw_gbm)

summary(logit_gbm)

logit_gbm_adjusted <- glm(rsv_outcome ~ covid + age_group + ce_week_days_numeric, 
                 data = covid_flu_iptw_results$dataset, family = quasibinomial(link="logit"), 
                 weights = iptw_gbm)

summary(logit_gbm_adjusted)

plot_summs(logit_adjustment_gbm, logit_no_weighted, logit_adjustment, logit_gbm_adjusted,
           model.names = c("GBM weights", "unweighted", "LR weights, adjusted", "GBM weights, adjusted"))

## For exporting the fit summary table comparison 
# export_summs(fit, fit2, scale = TRUE,
#              error_format = "[{conf.low}, {conf.high}]")

pairs_data <- covid_flu_iptw_results$dataset %>% 
  select(race_eth_cat, age_group, ce_week_days_numeric, sub_cohort) 

library(GGally)

ggpairs(pairs_data) +
  theme_bw()

## Need to make a logistic regression outcome model pipeline

### TODO have some way to output and save the logistic regression models 



#### Cox regression






