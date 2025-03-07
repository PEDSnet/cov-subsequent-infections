require(lubridate)
library(wesanderson)
library(survey)
library(survival)
library(WeightIt)
library(cobalt)
library(jtools)
require(cobalt)
require(WeightIt)
require(scattermore)
require(broom)
require(broom.mixed)
require(tableone)

# data RSV
message("sites with < 20% confirmed influenza test are excluded")
site_excl <- c("monte", "nationwide", "columbia", "emory", 
               "intermountain", "ochsner", "temple", "uth", 
               "osu", "utsw", "wakeforest")

postacute_min_start_date = as.Date("2022-04-01")
postacute_max_end_date = as.Date("2023-01-01")

cohort_entry_start_date = as.Date("2022-03-01")
cohort_entry_end_date = as.Date("2022-07-01")

cohort_1_label = "rsv_study"

comparison_group_string = "Influenza"


analytic_dataset <-
  results_tbl(paste0(cohort_1_label, "_analytic_highrisk")) %>% 
  rename(ce_date = ce_date.x) %>%
  filter(ce_date >= cohort_entry_start_date,
         ce_date < cohort_entry_end_date) %>%
  exclude_sites(site_excl)

analytic_dataset_final <-
  analytic_dataset %>% 
  filter(sex_cat != "Other/unknown") %>% 
  filter(is.na(covid_index_date_imputed) | covid_index_date_imputed==0) %>% 
  group_by(person_id) %>% 
  mutate(exclude = max(exclude_for_prior_rsv)) %>% 
  filter(exclude < 1) %>%
  slice_min(rsv_evidence_date, with_ties=FALSE) %>% 
  ungroup() %>% 
  mutate(lab_confirmed_index = ifelse(lab_confirmed==1, "lab confirmed", "dx only")) %>% 
  mutate(lab_confirmed_rsv_outcome = ifelse(lab_confirmed_rsv==1, "lab confirmed RSV", "dx only")) %>% 
  # mutate(visit_span_criteria_rank = as.character(visit_span_criteria_rank)) %>% 
  mutate(hospitalized_at_index = ifelse(hospital_flag==1 & visit_span_criteria_rank=="1", "Hospitalized", "Not hospitalized")) %>% 
  mutate(hospital_flag = as.character(hospital_flag)) %>% 
  mutate(hosp_with_vent = as.character(hosp_with_vent)) %>% 
  mutate(hospitalized_at_index = ifelse(is.na(hospitalized_at_index), "Not hospitalized", hospitalized_at_index)) %>% 
  mutate(hosp_with_vent = ifelse(is.na(hosp_with_vent), "0", hosp_with_vent))  %>% 
  mutate(hospitalized_at_index_detail = ifelse(hospitalized_at_index=="Hospitalized" & hosp_with_vent=="1", "Hospitalized with ventilation", hospitalized_at_index)) %>% collect() %>%
  mutate(pulmonary_respiratory = as.factor(pulmonary_respiratory),
         highrisk_flag = if_else(highrisk_flag, "Yes", "No"))

message("Figure 4a: RSV within 15-300 days (past post-acute)")
outcome_15_300 = "outcome_rsv_15_to_300"
comparison_group = "Influenza"
rsv_iptw_15_300_highrisk <- analytic_dataset_final %>%
  do_iptw(cohort_entry_start_date = cohort_entry_start_date,
          comparison_cohort = comparison_group,
          weight_formula = "covid~age_group+sex_cat+race_eth_cat+ce_week_factor+highrisk_flag",
          iptw_description = "covid_flu_outcome_rsv15_300_highrisk_site_excl")

df_15_300_highrisk <- results_tbl( paste0("covid_flu_outcome_rsv15_300_highrisk_site_excl",'_tabc')) %>%
  collect()
models_highrisk_4a <- list()
formula_no_pulmonary = "~ covid+ age_group + race_eth_cat + sex_cat + hospitalized_at_index +  util_inpatient +  util_outpatient +  util_ed + util_other"
models_highrisk_4a$unw_unadj_15_to_300 <- run_logistic_regression(df = df_15_300_highrisk,
                                                                  model_formula = paste0(outcome_15_300, " ~ covid"),
                                                                  weight_method = "none")

models_highrisk_4a$w_adj_15_to_300 <- run_logistic_regression(df = df_15_300_highrisk,
                                                              model_formula = paste0(outcome_15_300,
                                                                                     formula_no_pulmonary, "+ highrisk_flag"),
                                                              weight_method = "lr")
models_highrisk_4a$unw_adj_15_to_300 <- run_logistic_regression(df = df_15_300_highrisk,
                                                            model_formula = paste0(outcome_15_300,
                                                                                          formula_no_pulmonary, " + highrisk_flag"),
                                                            weight_method = "none")
models_highrisk_4a$w_unadj_15_to_300 <- run_logistic_regression(df = df_15_300_highrisk,
                                                      model_formula = paste0(outcome_15_300, " ~ covid"),
                                                      weight_method = "lr_s")
var_to_omit <- c("(Intercept)", "util_inpatientno_visits", "util_inpatientlow_utilizer", "util_inpatienthigh_utilizer",
                 "util_inpatientmoderate_utilizer",
                 "util_outpatientno_visits", "util_outpatientlow_utilizer", "util_outpatienthigh_utilizer", "util_outpatientmoderate_utilizer",
                 "util_edno_visits", "util_edlow_utilizer", "util_edhigh_utilizer", "util_edmoderate_utilizer",
                 "util_otherno_visits", "util_otherlow_utilizer", "util_otherhigh_utilizer", "util_othermoderate_utilizer")

p_4a <- plot_summs(models_highrisk_4a$unw_unadj_15_to_300,
                   models_highrisk_4a$unw_adj_15_to_300,
                   models_highrisk_4a$w_unadj_15_to_300,
                   models_highrisk_4a$w_adj_15_to_300,
                   omit.coefs = var_to_omit,
                   # models$model_weighted_adjusted_time,
                   model.names = c("Unweighted, unadjusted",
                                   "Unweighted, adjusted",
                                   "Weighted (LR), unadjusted",
                                   "Weighted (LR), adjusted"),
                   # "Weighted, adjusted for CE week"),
                   exp = TRUE)
p_4a <- p_4a + scale_y_discrete(labels = c("High Risk Conditions",
                                           "Hopspitalized for index infection",
                                           "Male",
                                           "Other/Unknown Race", "Non-Hispanic White", "Non-Hispanic Multiple", "Non-Hispanic Black/AA", "Non-Hispanic Asian/PI",
                                           "6m-1y old", "4-5y old", "3-4y old", "2-3y old", "1-2y old", "SARS-CoV-2 Infection")) +
  labs(x = "Odd ratios", y = "")
print(p_4a)
ggsave("results/figures_10082024/Figure_4a_rsv_15_300_highrisk_site_excl.png", plot = p_4a, width = 20, height = 10, units = "cm")





message("Figure 4b: RSV within 30-180 days (past post-acute)")
outcome_30_180 = "outcome_rsv_30_to_180"
comparison_group = "Influenza"
rsv_iptw_30_180_highrisk <- analytic_dataset_final %>%
  do_iptw(cohort_entry_start_date = cohort_entry_start_date,
          comparison_cohort = comparison_group,
          weight_formula = "covid~age_group+sex_cat+race_eth_cat+ce_week_factor+highrisk_flag",
          iptw_description = "covid_flu_outcome_rsv30_180_highrisk_site_excl")

df_30_180_highrisk <- results_tbl( paste0("covid_flu_outcome_rsv30_180_highrisk_site_excl",'_tabc')) %>%
  collect()
models_highrisk_4b <- list()
formula_no_pulmonary = "~ covid+ age_group + race_eth_cat + sex_cat + hospitalized_at_index +  util_inpatient +  util_outpatient +  util_ed + util_other"
models_highrisk_4b$unw_unadj_30_to_180 <- run_logistic_regression(df = df_30_180_highrisk,
                                                                  model_formula = paste0(outcome_30_180, " ~ covid"),
                                                                  weight_method = "none")

models_highrisk_4b$w_adj_30_to_180 <- run_logistic_regression(df = df_30_180_highrisk,
                                                              model_formula = paste0(outcome_30_180,
                                                                                     formula_no_pulmonary, "+ highrisk_flag"),
                                                              weight_method = "lr")
models_highrisk_4b$unw_adj_30_to_180 <- run_logistic_regression(df = df_30_180_highrisk,
                                                                model_formula = paste0(outcome_30_180,
                                                                                       formula_no_pulmonary, " + highrisk_flag"),
                                                                weight_method = "none")
models_highrisk_4b$w_unadj_30_to_180 <- run_logistic_regression(df = df_30_180_highrisk,
                                                                model_formula = paste0(outcome_30_180, " ~ covid"),
                                                                weight_method = "lr_s")
var_to_omit <- c("(Intercept)", "util_inpatientno_visits", "util_inpatientlow_utilizer", "util_inpatienthigh_utilizer",
                 "util_inpatientmoderate_utilizer",
                 "util_outpatientno_visits", "util_outpatientlow_utilizer", "util_outpatienthigh_utilizer", "util_outpatientmoderate_utilizer",
                 "util_edno_visits", "util_edlow_utilizer", "util_edhigh_utilizer", "util_edmoderate_utilizer",
                 "util_otherno_visits", "util_otherlow_utilizer", "util_otherhigh_utilizer", "util_othermoderate_utilizer")

p_4b <- plot_summs(models_highrisk_4b$unw_unadj_30_to_180,
                   models_highrisk_4b$unw_adj_30_to_180,
                   models_highrisk_4b$w_unadj_30_to_180,
                   models_highrisk_4b$w_adj_30_to_180,
                   omit.coefs = var_to_omit,
                   # models$model_weighted_adjusted_time,
                   model.names = c("Unweighted, unadjusted",
                                   "Unweighted, adjusted",
                                   "Weighted (LR), unadjusted",
                                   "Weighted (LR), adjusted"),
                   # "Weighted, adjusted for CE week"),
                   exp = TRUE)
p_4b <- p_4b + scale_y_discrete(labels = c("High Risk Conditions",
                                           "Hopspitalized for index infection",
                                           "Male",
                                           "Other/Unknown Race", "Non-Hispanic White", "Non-Hispanic Multiple", "Non-Hispanic Black/AA", "Non-Hispanic Asian/PI",
                                           "6m-1y old", "4-5y old", "3-4y old", "2-3y old", "1-2y old", "SARS-CoV-2 Infection")) +
  labs(x = "Odd ratios", y = "")
print(p_4b)
ggsave("results/figures_10082024/Figure_4b_rsv_30_180_highrisk_site_excl.png", plot = p_4b, width = 20, height = 10, units = "cm")

p_s3a <- rsv_iptw_30_180_highrisk$love_plot
ggsave("results/figures_10082024/Figure_s3a_rsv_30_180_highrisk_site_excl.png",p_s3a,  width = 6.5, height = 10, unit = "in")

message("Figure 4a: RSV within 0-180 days (include co-infections)")
# create this outcome from 0_to_60 and 61_to_180 outcomes
analytic_dataset_final <- analytic_dataset_final%>%
                              mutate(outcome_rsv_0_to_180 = ifelse(outcome_rsv_0_to_60 == 1 | outcome_rsv_61_to_180 == 1, 1, 0))

outcome_0_180 = "outcome_rsv_0_to_180"
comparison_group = "Influenza"
rsv_iptw_0_180_highrisk <- analytic_dataset_final %>%
  do_iptw(cohort_entry_start_date = cohort_entry_start_date,
          comparison_cohort = comparison_group,
          weight_formula = "covid~age_group+sex_cat+race_eth_cat+ce_week_factor+highrisk_flag",
          iptw_description = "covid_flu_outcome_rsv0_180_highrisk_site_excl")

df_0_180_highrisk <- results_tbl( paste0("covid_flu_outcome_rsv0_180_highrisk_site_excl",'_tabc')) %>%
  collect()
models_highrisk_4a <- list()
formula_no_pulmonary = "~ covid+ age_group + race_eth_cat + sex_cat + hospitalized_at_index +  util_inpatient +  util_outpatient +  util_ed + util_other"
models_highrisk_4a$unw_unadj_0_to_180 <- run_logistic_regression(df = df_0_180_highrisk,
                                                                  model_formula = paste0(outcome_0_180, " ~ covid"),
                                                                  weight_method = "none")

models_highrisk_4a$w_adj_0_to_180 <- run_logistic_regression(df = df_0_180_highrisk,
                                                              model_formula = paste0(outcome_0_180,
                                                                                     formula_no_pulmonary, "+ highrisk_flag"),
                                                              weight_method = "lr")
models_highrisk_4a$unw_adj_0_to_180 <- run_logistic_regression(df = df_0_180_highrisk,
                                                                model_formula = paste0(outcome_0_180,
                                                                                       formula_no_pulmonary, " + highrisk_flag"),
                                                                weight_method = "none")
models_highrisk_4a$w_unadj_0_to_180 <- run_logistic_regression(df = df_0_180_highrisk,
                                                                model_formula = paste0(outcome_0_180, " ~ covid"),
                                                                weight_method = "lr_s")
var_to_omit <- c("(Intercept)", "util_inpatientno_visits", "util_inpatientlow_utilizer", "util_inpatienthigh_utilizer",
                 "util_inpatientmoderate_utilizer",
                 "util_outpatientno_visits", "util_outpatientlow_utilizer", "util_outpatienthigh_utilizer", "util_outpatientmoderate_utilizer",
                 "util_edno_visits", "util_edlow_utilizer", "util_edhigh_utilizer", "util_edmoderate_utilizer",
                 "util_otherno_visits", "util_otherlow_utilizer", "util_otherhigh_utilizer", "util_othermoderate_utilizer")

p_4a <- plot_summs(models_highrisk_4a$unw_unadj_0_to_180,
                   models_highrisk_4a$unw_adj_0_to_180,
                   models_highrisk_4a$w_unadj_0_to_180,
                   models_highrisk_4a$w_adj_0_to_180,
                   omit.coefs = var_to_omit,
                   # models$model_weighted_adjusted_time,
                   model.names = c("Unweighted, unadjusted",
                                   "Unweighted, adjusted",
                                   "Weighted (LR), unadjusted",
                                   "Weighted (LR), adjusted"),
                   # "Weighted, adjusted for CE week"),
                   exp = TRUE)
p_4a <- p_4a + scale_y_discrete(labels = c("High Risk Conditions",
                                           "Hopspitalized for index infection",
                                           "Male",
                                           "Other/Unknown Race", "Non-Hispanic White", "Non-Hispanic Multiple", "Non-Hispanic Black/AA", "Non-Hispanic Asian/PI",
                                           "6m-1y old", "4-5y old", "3-4y old", "2-3y old", "1-2y old", "SARS-CoV-2 Infection")) +
  labs(x = "Odd ratios", y = "")
print(p_4a)
ggsave("results/figures_10082024/Figure_4a_rsv_0_180_highrisk_site_excl.png", plot = p_4a, width = 20, height = 10, units = "cm")

# p_s3a <- rsv_iptw_30_180_highrisk$love_plot
# ggsave("results/figures_10082024/Figure_s3a_rsv_30_180_highrisk_site_excl.png",p_s3a,  width = 6.5, height = 10, unit = "in")
