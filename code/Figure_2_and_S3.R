# this script generates Figure 2 and S3 for the manuscript

require(tibble)
require(knitr)
require(kableExtra)
require(readr)
require(table1)
require(ggplot2)
library(GGally)
require(tidyverse)
require(lubridate)
library(wesanderson)
library(survey)
library(survival)
library(WeightIt)
library(cobalt)
library(jtools)
require(cobalt)
require(WeightIt)
# require(scattermore)
require(broom)
require(broom.mixed)
require(tableone)
require(forestplot)
# Data RSV
postacute_min_start_date = as.Date("2022-04-01")
postacute_max_end_date = as.Date("2023-01-01")

cohort_entry_start_date = as.Date("2022-03-01")
cohort_entry_end_date = as.Date("2022-07-01")

cohort_1_label = "rsv_study"

comparison_group_string = "Influenza"

message("sites with < 20% confirmed influenza test are excluded")
site_excl <- c("monte", "nationwide", "columbia", "emory", 
               "intermountain", "ochsner", "temple", "uth", 
               "osu", "utsw", "wakeforest")

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

# Data Resp
postacute_min_start_date_resp = as.Date("2021-07-01")
postacute_max_end_date_resp = as.Date("2023-03-01")

cohort_entry_start_date_resp = as.Date("2021-06-01")
cohort_entry_end_date_resp = as.Date("2023-01-01")
cohort_entry_end_date_resp = as.Date("2022-10-01")

cohort_1_label_resp = "resp_study"


# description = "covid_flu_3"
# description = "covid_flu_noutil"
# description = "covid_flu_util"
# description = "covid_ari_util"
# description = "covid_ari_ce_week"
# description = "covid_flu_nonexclude"

description_resp = "covid_flu_comparison"
# description = "covid_ari_comparison"
# comparison_group_string = "Respiratory"
comparison_group_string_resp = "Influenza"


analytic_dataset_resp <-
  results_tbl(paste0(cohort_1_label_resp, "_analytic_highrisk")) %>% 
  filter(ce_date >= cohort_entry_start_date_resp,
         ce_date < cohort_entry_end_date_resp) %>%
  exclude_sites(site_excl)

analytic_dataset_final_resp <-
  analytic_dataset_resp %>% 
  filter(sex_cat != "Other/unknown") %>% 
  filter(is.na(covid_index_date_imputed) | covid_index_date_imputed==0) %>% 
  group_by(person_id) %>% # TODO might need to filter earliest resp date or something
  ungroup() %>% 
  mutate(lab_confirmed_index = ifelse(lab_confirmed==1, "lab confirmed", "dx only")) %>% 
  # mutate(visit_span_criteria_rank = as.character(visit_span_criteria_rank)) %>% 
  mutate(hospitalized_at_index = if_else(hospital_flag==1 & visit_span_criteria_rank=="1", "Yes","No")) %>% 
  mutate(hospital_flag = as.character(hospital_flag)) %>% 
  mutate(hosp_with_vent = as.character(hosp_with_vent)) %>% 
  mutate(hospitalized_at_index = if_else(is.na(hospitalized_at_index), "No", hospitalized_at_index)) %>% 
  mutate(hosp_with_vent = ifelse(is.na(hosp_with_vent), "0", hosp_with_vent))  %>% 
  mutate(hospitalized_at_index_detail = ifelse(hospitalized_at_index== "Yes" & hosp_with_vent=="1", "Hospitalized with ventilation", hospitalized_at_index)) %>% 
  mutate(resp_outcome_is_reinfection = case_when(resp_outcome_flu_related=="x" & sub_cohort=="Influenza" ~ "1",
                                                 outcome_covid30d=="1" & sub_cohort=="COVID" ~ "1",
                                                 TRUE ~ "0")) %>% 
  mutate(highrisk_flag = if_else(highrisk_flag, "Yes", "No")) %>%
  # mutate(rsv_outcome_postacute = ifelse(rsv_occurrence_period == "post acute (15-180 days after index)", 1, 0)) %>% 
  collect()


# New Figure 2 (Old F7): Sensitivity analysis RSV outcome, high risk conditions
## New Fig 2a: RSV 3 periods, 4 models 
models_highrisk_2a = list()

outcome_0_60 = "outcome_rsv_0_to_60"
comparison_group = "Influenza"

# formula_no_pulmonary = "~ `SARS-CoV-2 Infection: ` + `Age groups: ` + `Race: ` + `Gender: ` + `Hospitalized for index infection: ` +  util_inpatient +  util_outpatient +  util_ed + util_other"
formula_no_pulmonary = "~ covid+ age_group + race_eth_cat + sex_cat + hospitalized_at_index +  util_inpatient +  util_outpatient +  util_ed + util_other"
iptw_0_60_highrisk <- analytic_dataset_final %>%
  do_iptw(cohort_entry_start_date = cohort_entry_start_date,
          comparison_cohort = comparison_group,
          weight_formula = "covid~age_group+sex_cat+race_eth_cat+ce_week_factor+highrisk_flag",
          iptw_description = "covid_flu_outcome_rsv0_60_highrisk_site_ex")
df_0_60_highrisk <- results_tbl( paste0("covid_flu_outcome_rsv0_60_highrisk_site_ex", "_tabc")) %>% 
  collect()

models_highrisk_2a$unw_unadj_0_to_60 <- run_logistic_regression(df = df_0_60_highrisk,
                                                                model_formula = paste0(outcome_0_60, " ~ covid"),
                                                                weight_method = "none")

models_highrisk_2a$w_adj_0_to_60 <- run_logistic_regression(df = df_0_60_highrisk,
                                                            model_formula = paste0(outcome_0_60,
                                                                                   formula_no_pulmonary, " + highrisk_flag"), 
                                                            weight_method = "lr")
outcome_1_60 = "outcome_rsv_1_to_60"
iptw_1_60_highrisk <- analytic_dataset_final %>%
  do_iptw(cohort_entry_start_date = cohort_entry_start_date,
          comparison_cohort = comparison_group,
          weight_formula = "covid~age_group+sex_cat+race_eth_cat+ce_week_factor+highrisk_flag",
          iptw_description = "covid_flu_outcome_rsv1_60_highrisk_site_ex")

df_1_60_highrisk <- results_tbl( paste0("covid_flu_outcome_rsv1_60_highrisk_site_ex",'_tabc')) %>% 
  collect()

models_highrisk_2a$unw_unadj_1_to_60 <- run_logistic_regression(df = df_1_60_highrisk,
                                                                model_formula = paste0(outcome_1_60, " ~ covid"),
                                                                weight_method = "none")

models_highrisk_2a$w_adj_1_to_60 <- run_logistic_regression(df = df_1_60_highrisk,
                                                            model_formula = paste0(outcome_1_60,
                                                                                   formula_no_pulmonary, " + highrisk_flag"), 
                                                            weight_method = "lr")

outcome_15_180 = "outcome_rsv_15_to_180"

rsv_iptw_15_180_highrisk <- analytic_dataset_final %>%
  do_iptw(cohort_entry_start_date = cohort_entry_start_date,
          comparison_cohort = comparison_group,
          weight_formula = "covid~age_group+sex_cat+race_eth_cat+ce_week_factor+highrisk_flag",
          iptw_description = "covid_flu_outcome_rsv15_180_highrisk_site_excl")

df_15_180_highrisk <- results_tbl( paste0("covid_flu_outcome_rsv15_180_highrisk_site_excl",'_tabc')) %>%
  collect()

models_highrisk_2a$unw_unadj_15_to_180 <- run_logistic_regression(df = df_15_180_highrisk,
                                                                  model_formula = paste0(outcome_15_180, " ~ covid"),
                                                                  weight_method = "none")

models_highrisk_2a$w_adj_15_to_180 <- run_logistic_regression(df = df_15_180_highrisk,
                                                              model_formula = paste0(outcome_15_180,
                                                                                     formula_no_pulmonary, "+ highrisk_flag"),
                                                              weight_method = "lr")

var_to_omit <- c("(Intercept)", "util_inpatientno_visits", "util_inpatientlow_utilizer", "util_inpatienthigh_utilizer",
                 "util_inpatientmoderate_utilizer",
                 "util_outpatientno_visits", "util_outpatientlow_utilizer", "util_outpatienthigh_utilizer", "util_outpatientmoderate_utilizer",
                 "util_edno_visits", "util_edlow_utilizer", "util_edhigh_utilizer", "util_edmoderate_utilizer",
                 "util_otherno_visits", "util_otherlow_utilizer", "util_otherhigh_utilizer", "util_othermoderate_utilizer")

p_2a_3periods <- plot_summs(models_highrisk_2a$unw_unadj_0_to_60,
                            models_highrisk_2a$w_adj_0_to_60,
                            models_highrisk_2a$unw_unadj_1_to_60,
                            models_highrisk_2a$w_adj_1_to_60,
                            models_highrisk_2a$unw_unadj_15_to_180,
                            models_highrisk_2a$w_adj_15_to_180,
                            omit.coefs = var_to_omit,
                            model.names = c("RSV 0-60d: Unweighted, unadjusted",
                                            "RSV 0-60d: Weighted, adjusted",
                                            "RSV 1-60d: Unweighted, unadjusted",
                                            "RSV 1-60d: Weighted, adjusted",
                                            "RSV 15-180d: Unweighted, unadjusted",
                                            "RSV 15-180d: Weighted, adjusted"),
                            # "Weighted, adjusted for CE week"),
                            exp = TRUE)

p_2a_3periods <- p_2a_3periods + scale_y_discrete(labels = c("High risk conditions",
                                                             "Hopspitalized for index infection",
                                                             "Male",
                                                             "Other/Unknown Race", "Non-Hispanic White", "Non-Hispanic Multiple", "Non-Hispanic Black/AA", "Non-Hispanic Asian/PI",
                                                             "6m-1y old", "4-5y old", "3-4y old", "2-3y old", "1-2y old", "SARS-CoV-2 Infection")) +
  labs(x = "Odd ratios", y = "")
print(p_2a_3periods)
ggsave("results/figures_10082024/Figure_2a_rsv_highrisk_3periods_site_excl.png", plot = p_2a_3periods, width = 20, height = 10, unit = "cm")


## New Fig 2a: RSV 1 period, 4 models 
models_highrisk_2a$unw_adj_15_to_180 <- run_logistic_regression(df = df_15_180_highrisk,
                                                                model_formula = paste0(outcome_15_180, formula_no_pulmonary, "+ highrisk_flag"),
                                                                weight_method = "none")

models_highrisk_2a$w_unadj_15_to_180 <- run_logistic_regression(df = df_15_180_highrisk,
                                                                model_formula = paste0(outcome_15_180,
                                                                                       "~ covid"),
                                                                weight_method = "lr")

# var_to_omit <- c("(Intercept)", "util_inpatientno_visits", "util_inpatientlow_utilizer", "util_inpatienthigh_utilizer",
#                  "util_inpatientmoderate_utilizer",
#                  "util_outpatientno_visits", "util_outpatientlow_utilizer", "util_outpatienthigh_utilizer", "util_outpatientmoderate_utilizer",
#                  "util_edno_visits", "util_edlow_utilizer", "util_edhigh_utilizer", "util_edmoderate_utilizer",
#                  "util_otherno_visits", "util_otherlow_utilizer", "util_otherhigh_utilizer", "util_othermoderate_utilizer")
# 
p_2a <- plot_summs(models_highrisk_2a$unw_unadj_15_to_180,
                   models_highrisk_2a$unw_adj_15_to_180,
                   models_highrisk_2a$w_unadj_15_to_180,
                   models_highrisk_2a$w_adj_15_to_180,
                   omit.coefs = var_to_omit,
                   model.names = c("Unweighted, unadjusted",
                                   "Unweighted, adjusted",
                                   "Weighted (LR), unadjusted",
                                   "Weighted (LR), adjusted"),
                   # "Weighted, adjusted for CE week"),
                   exp = TRUE)

p_2a <- p_2a + scale_y_discrete(labels = c("High Risk Conditions",
                                           "Hopspitalized for index infection",
                                           "Male",
                                           "Other/Unknown Race", "Non-Hispanic White", "Non-Hispanic Multiple", "Non-Hispanic Black/AA", "Non-Hispanic Asian/PI",
                                           "6m-1y old", "4-5y old", "3-4y old", "2-3y old", "1-2y old", "SARS-CoV-2 Infection")) +
  labs(x = "Odd ratios", y = "")
print(p_2a)
ggsave("results/figures_10082024/Figure_2a_rsv_15_180_highrisk_1period_site_excl.png", plot = p_2a, width = 20, height = 10, unit = "cm")

export_summs(models_highrisk_2a$unw_unadj_15_to_180,
             models_highrisk_2a$unw_adj_15_to_180,
             models_highrisk_2a$w_unadj_15_to_180,
             models_highrisk_2a$w_adj_15_to_180,
             omit.coefs = var_to_omit,
             model.names = c("Unweighted, unadjusted",
                             "Unweighted, adjusted",
                             "Weighted (LR), unadjusted",
                             "Weighted (LR), adjusted"),
             # "Weighted, adjusted for CE week"),
             exp = TRUE,
             error_format = "[{conf.low}, {conf.high}]",
             to.file = "html", file.name = "fig2a_1period_stats.html")


# Fig 2b: Subsequent Respiratory infection within 30-180 days **Including covid, or flu, ANY respiratory infection within 30-180 days**
models_highrisk_2b = list()

outcome_ari = "resp30d_or_covid30d"
comparison_group_ari = "Influenza"

formula_no_pulmonary = "~ covid+ age_group + race_eth_cat + sex_cat + hospitalized_at_index +  util_inpatient +  util_outpatient +  util_ed + util_other"

ari_iptw_30_180 <- analytic_dataset_final_resp %>%
  mutate(ce_month = floor_date(ce_date, "month")) %>%
  mutate(ce_month_s = as.character(ce_month)) %>%
  filter(resp_outcome_is_reinfection=="0") %>% # TODO allow for re-infection
  mutate(resp30d_or_covid30d = case_when(outcome_resp30d==1 ~ 1,
                                         outcome_covid30d==1 ~ 1,
                                         TRUE ~ 0)) %>%
  mutate(covid = ifelse(sub_cohort == "COVID", 1, 0)) %>%
  do_iptw(cohort_entry_start_date = cohort_entry_start_date,
          comparison_cohort = comparison_group,
          weight_formula = "covid~age_group+sex_cat+race_eth_cat+pulmonary_respiratory+hospitalized_at_index+ce_week_factor+ highrisk_flag", ## PULMONARY_RESP?
          iptw_description = "resp_notime_w_flu_site_exc_resp30d_highrisk")
dataset_iptw_ari <- results_tbl( paste0("resp_notime_w_flu_site_exc_resp30d_highrisk",'_tabc')) %>%
  collect()

models_ari_2b = list()

models_ari_2b$model_unweighted_unadjusted <- run_logistic_regression(df = dataset_iptw_ari,
                                                                     model_formula = paste0(outcome_ari, " ~ covid"),
                                                                     weight_method = "none")

models_ari_2b$model_unweighted_adjusted <- run_logistic_regression(df = dataset_iptw_ari,
                                                                   model_formula = paste0(outcome_ari,
                                                                                          formula_no_pulmonary, " + highrisk_flag"),
                                                                   weight_method = "none")

models_ari_2b$model_weighted_unadjusted <- run_logistic_regression(df = dataset_iptw_ari,
                                                                   model_formula = paste0(outcome_ari, " ~ covid"),
                                                                   weight_method = "lr_s")

# # TODO something is wrong with inpatient
models_ari_2b$model_weighted_adjusted <- run_logistic_regression(df = dataset_iptw_ari,
                                                                 model_formula = paste0(outcome_ari,
                                                                                        formula_no_pulmonary, " + highrisk_flag"),
                                                                 weight_method = "lr_s")

var_to_omit <- c("(Intercept)", "util_inpatientno_visits", "util_inpatientlow_utilizer", "util_inpatienthigh_utilizer",
                 "util_inpatientmoderate_utilizer",
                 "util_outpatientno_visits", "util_outpatientlow_utilizer", "util_outpatienthigh_utilizer", "util_outpatientmoderate_utilizer",
                 "util_edno_visits", "util_edlow_utilizer", "util_edhigh_utilizer", "util_edmoderate_utilizer",
                 "util_otherno_visits", "util_otherlow_utilizer", "util_otherhigh_utilizer", "util_othermoderate_utilizer")

p_2b <- plot_summs(models_ari_2b$model_unweighted_unadjusted,
                   models_ari_2b$model_unweighted_adjusted,
                   models_ari_2b$model_weighted_unadjusted,
                   models_ari_2b$model_weighted_adjusted,
                   omit.coefs = var_to_omit,
                   # models$model_weighted_adjusted_time,
                   model.names = c("Unweighted, unadjusted",
                                   "Unweighted, adjusted",
                                   "Weighted (LR), unadjusted",
                                   "Weighted (LR), adjusted"),
                   # "Weighted, adjusted for CE week"),
                   exp = TRUE)
p_2b <- p_2b + scale_y_discrete(labels = c("High Risk Conditions",
                                           "Hopspitalized for index infection",
                                           "Male",
                                           "Other/Unknown Race", "Non-Hispanic White", "Non-Hispanic Multiple", "Non-Hispanic Black/AA", "Non-Hispanic Asian/PI",
                                           "6m-1y old", "4-5y old", "3-4y old", "2-3y old", "1-2y old", "SARS-CoV-2 Infection")) +
  labs(x = "Odd ratios", y = "")
print(p_2b)
ggsave("results/figures_10082024/Figure_2b_ari_with_covid_30_180_highrisk_site_excl.png", plot = p_2b, width = 20, height = 10, units = "cm")

export_summs(models_ari_2b$model_unweighted_unadjusted,
             models_ari_2b$model_unweighted_adjusted,
             models_ari_2b$model_weighted_unadjusted,
             models_ari_2b$model_weighted_adjusted,
             # models$model_weighted_adjusted_time,
             model.names = c("Unweighted, unadjusted",
                             "Unweighted, adjusted",
                             "Weighted (LR), unadjusted",
                             "Weighted (LR), adjusted"),
             # "Weighted, adjusted for CE week"),
             exp = TRUE,
             error_format = "[{conf.low}, {conf.high}]",
             to.file = "html", file.name = "fig_2b_stats.html")


models_highrisk_2c = list()

outcome = "outcome_general30d"
comparison_group = "Influenza"

any_inf_iptw_flu_30_180 <- analytic_dataset_final_resp %>%
  filter(resp_outcome_is_reinfection=="0") %>% # TODO allow for re-infection
  do_iptw(cohort_entry_start_date = cohort_entry_start_date_resp,
          comparison_cohort = comparison_group,
          weight_formula = "covid~age_group+sex_cat+race_eth_cat+ce_week_factor+highrisk_flag",
          iptw_description = "resp_study_covid_flu_site_exc_general30d_highrisk")

dataset_with_iptw <- results_tbl( paste0("resp_study_covid_flu_site_exc_general30d_highrisk", '_tabc')) %>%
  collect()

models_highrisk_2c$model_unweighted_unadjusted <- run_logistic_regression(df = dataset_with_iptw,
                                                                          model_formula = paste0(outcome, " ~ covid"),
                                                                          weight_method = "none")

models_highrisk_2c$model_unweighted_adjusted <- run_logistic_regression(df = dataset_with_iptw,
                                                                        model_formula = paste0(outcome,
                                                                                               formula_no_pulmonary, " + highrisk_flag"), 
                                                                        weight_method = "none")

models_highrisk_2c$model_weighted_unadjusted <- run_logistic_regression(df = dataset_with_iptw,
                                                                        model_formula = paste0(outcome, " ~ covid"),
                                                                        weight_method = "lr")

models_highrisk_2c$model_weighted_adjusted <- run_logistic_regression(df = dataset_with_iptw,
                                                                      model_formula = paste0(outcome,
                                                                                             formula_no_pulmonary, " + highrisk_flag"), 
                                                                      weight_method = "lr")
var_to_omit <- c("(Intercept)", "util_inpatientno_visits", "util_inpatientlow_utilizer", "util_inpatienthigh_utilizer",
                 "util_inpatientmoderate_utilizer",
                 "util_outpatientno_visits", "util_outpatientlow_utilizer", "util_outpatienthigh_utilizer", "util_outpatientmoderate_utilizer",
                 "util_edno_visits", "util_edlow_utilizer", "util_edhigh_utilizer", "util_edmoderate_utilizer",
                 "util_otherno_visits", "util_otherlow_utilizer", "util_otherhigh_utilizer", "util_othermoderate_utilizer")

p_2c <- plot_summs(models_highrisk_2c$model_unweighted_unadjusted,
                   models_highrisk_2c$model_unweighted_adjusted,
                   models_highrisk_2c$model_weighted_unadjusted,
                   models_highrisk_2c$model_weighted_adjusted,
                   omit.coefs = var_to_omit,
                   # models$model_weighted_adjusted_time,
                   model.names = c("Unweighted, unadjusted",
                                   "Unweighted, adjusted",
                                   "Weighted (LR), unadjusted",
                                   "Weighted (LR), adjusted"),
                   # "Weighted, adjusted for CE week"),
                   exp = TRUE)

p_2c <- p_2c + scale_y_discrete(labels = c("High risk conditions", 
                                           "Hopspitalized for index infection", 
                                           "Male", 
                                           "Other/Unknown Race", "Non-Hispanic White", "Non-Hispanic Multiple", "Non-Hispanic Black/AA", "Non-Hispanic Asian/PI",
                                           "6m-1y old", "4-5y old", "3-4y old", "2-3y old", "1-2y old", "SARS-CoV-2 Infection")) +
  labs(x = "Odd ratios", y = "") 

ggsave("results/figures_10082024/Figure_2c_general_30_180_highrisk_site_excl.png", plot = p_2c, width = 20, height = 10, units = "cm")

# export_summs(models_highrisk_2c$model_unweighted_unadjusted,
#            models_highrisk_2c$model_unweighted_adjusted,
#            models_highrisk_2c$model_weighted_unadjusted,
#            models_highrisk_2c$model_weighted_adjusted,
#            # models$model_weighted_adjusted_time,
#            model.names = c("Unweighted, unadjusted",
#                            "Unweighted, adjusted",
#                            "Weighted (LR), unadjusted",
#                            "Weighted (LR), adjusted"),
#                            # "Weighted, adjusted for CE week"),
#            exp = TRUE,
#            error_format = "[{conf.low}, {conf.high}]")


# Figure S3
p_s3a <- rsv_iptw_15_180_highrisk$love_plot
ggsave("results/figures_10082024/Figure_s3a_rsv_15_180_highrisk_site_excl.png",p_s3a,  width = 6.5, height = 10, unit = "in")

p_s3b <- ari_iptw_30_180$love_plot
ggsave("results/figures_10082024/Figure_s3b_ari_30_180_highrisk_site_excl.png",p_s3b,  width = 6.5, height = 10, unit = "in")

p_s3c <- any_inf_iptw_flu_30_180$love_plot
ggsave("results/figures_10082024/Figure_s3c_general_30_180_highrisk_site_excl.png", p_s3c, width = 6.5, height = 10, unit = "in")


## New Fig 2a: RSV 15-180, 1 model: weighted adjusted
outcome_15_180 = "outcome_rsv_15_to_180"
comparison_group = "Influenza"

formula_no_pulmonary = "~ covid+ age_group + race_eth_cat + sex_cat + hospitalized_at_index +  util_inpatient +  util_outpatient +  util_ed + util_other"

rsv_iptw_15_180_highrisk <- analytic_dataset_final %>%
  do_iptw(cohort_entry_start_date = cohort_entry_start_date,
          comparison_cohort = comparison_group,
          weight_formula = "covid~age_group+sex_cat+race_eth_cat+ce_week_factor",
          iptw_description = "covid_flu_rsv15_180_highrisk_site_excl")

df_15_180_highrisk <- results_tbl( paste0("covid_flu_rsv15_180_highrisk_site_excl",'_tabc')) %>%
  collect() %>% 
  mutate(hospitalized_at_index = factor(hospitalized_at_index, levels = c("Not hospitalized", "Hospitalized"))) %>%
  mutate(age_group = factor(age_group, levels = c("<6 months", "6m-1y old", "1-2y old", "2-3y old", "3-4y old", "4-5y old")))

# weighted, adjusted models
models_highrisk_2a <- run_logistic_regression(df = df_15_180_highrisk,
                                              model_formula = paste0(outcome_15_180, formula_no_pulmonary, "+ highrisk_flag"),
                                              weight_method = "lr")

base_data_2a <- as.data.frame(coef(models_highrisk_2a))
base_data_2a <- cbind.data.frame(base_data_2a, confint(models_highrisk_2a))
base_data_2a <- cbind.data.frame(base_data_2a, as.data.frame(coef(summary(models_highrisk_2a))[,4]))

var_to_omit <- c("(Intercept)", "util_inpatientno_visits", "util_inpatientlow_utilizer", "util_inpatienthigh_utilizer",
                 "util_inpatientmoderate_utilizer",
                 "util_outpatientno_visits", "util_outpatientlow_utilizer", "util_outpatienthigh_utilizer", "util_outpatientmoderate_utilizer",
                 "util_edno_visits", "util_edlow_utilizer", "util_edhigh_utilizer", "util_edmoderate_utilizer",
                 "util_otherno_visits", "util_otherlow_utilizer", "util_otherhigh_utilizer", "util_othermoderate_utilizer")

base_data_2a <- base_data_2a %>% rename(OR = "coef(models_highrisk_2a)", 
                                   pval = "coef(summary(models_highrisk_2a))[, 4]", 
                                   lower = "2.5 %", upper = "97.5 %") %>%
  mutate(covariates = rownames(.)) %>%
  filter(!(covariates %in% var_to_omit)) %>%
  mutate(across(c("OR", "lower", "upper"), ~exp(.x))) %>%
  mutate(across(c("OR", "lower", "upper"), ~round(.x, 2))) %>%
  mutate(ci = paste0(lower, "-", upper)) %>%
  mutate(covariates = case_when(covariates == "highrisk_flagYes"~ "At least 1 high risk condition", 
                                covariates == "hospitalized_at_indexYes" ~ "Hopspitalized for index infection",
                                covariates == "race_eth_catOther/Unknown" ~ "Other/Unknown race",
                                covariates == "race_eth_catHispanic" ~ "Hispanic",
                                covariates == "race_eth_catNH_White" ~ "Non-Hispanic White",
                                covariates == "race_eth_catNH_Black/AA" ~ "Non-Hispanic Black/AA",
                                covariates == "race_eth_catNH_Asian/PI" ~ "Non-Hispanic Asian/PI",
                                covariates == "race_eth_catNH_Multiple" ~ "Non-Hispanic Multiple races",
                                covariates == "sex_catMale" ~ "Male",
                                covariates == "sex_catFemale" ~ "Female",
                                covariates == "age_group<6 months" ~ "<6 months",
                                covariates == "age_group6m-1y old" ~ "6m-1y old",
                                covariates == "age_group1-2y old" ~ "1-2y old",
                                covariates == "age_group2-3y old" ~ "2-3y old",
                                covariates == "age_group3-4y old" ~ "3-4y old",
                                covariates == "age_group4-5y old" ~ "4-5y old",
                                covariates == "covid" ~ "SARS-CoV-2 Infections",
                                TRUE ~ covariates))%>%
  mutate(pval = case_when(pval > 0.05 ~ format(round(pval, 4), scientific = FALSE),
                          pval >= 0.005 ~ paste0(format(round(pval, 4), scientific = FALSE), "*"),
                          pval >= 0.001 ~ paste0(format(round(pval, 4), scientific = FALSE), "**"),
                          pval >= 0.0001 ~ paste0(format(round(pval, 4), scientific = FALSE), "***"),
                          TRUE ~ "<0.0001****"))

png("results/figures_10082024/Figure_2a_rsv15_180_highrisk_site_excl_w_adj.png", width = 20, height = 10, units = "cm", res = 300 )
p2a_w_adj <- base_data_2a |>
  forestplot(labeltext = c(covariates, OR, ci, pval),
             mean = OR, 
             clip = c(0.1, 5.0), #xlim
             xlog = FALSE,
             xlab = "Odds ratios", boxsize = 0.25, zero = 1) |>
  fp_set_style(box = "darkblue",
               line = "darkblue") |>
  fp_add_header(covariates = c("", ""),
                ci = c("95% CI", ""),
                OR = c("OR", ""),
                pval = c("p-value", ""))      
print(p2a_w_adj)

dev.off()

p_s3a <- rsv_iptw_15_180_highrisk$love_plot + 
  geom_vline(xintercept = 0.2, linetype = "dashed") +
  geom_vline(xintercept = -0.2, linetype = "dashed")
print(p_s3a)
ggsave("results/figures_10082024/Figure_s3a_rsv_15_180_highrisk_site_excl.png",p_s3a,  width = 6.5, height = 10, unit = "in")


## New Fig 2b: ari any resp infections 15-180, 1 model: weighted adjusted
outcome_ari = "resp15d_or_covid15d"
comparison_group_ari = "Influenza"

formula_no_pulmonary = "~ covid+ age_group + race_eth_cat + sex_cat + hospitalized_at_index +  util_inpatient +  util_outpatient +  util_ed + util_other"

ari_iptw_15_180 <- analytic_dataset_final_resp %>%
  mutate(ce_month = floor_date(ce_date, "month")) %>%
  mutate(ce_month_s = as.character(ce_month)) %>%
  filter(resp_outcome_is_reinfection=="0") %>% # TODO allow for re-infection
  mutate(resp15d_or_covid15d = case_when(outcome_resp15d==1 ~ 1,
                                         outcome_covid15d==1 ~ 1,
                                         TRUE ~ 0)) %>%
  #mutate(hospitalized_at_index = ifelse(hospitalized_at_index == "Yes", 1, 0)) %>%
  #mutate(highrisk_flag = ifelse(highrisk_flag == "Yes", 1, 0))%>%
  do_iptw(#analytic_dataset = ari_iptw_15_180, 
    cohort_entry_start_date = cohort_entry_start_date_resp,
    comparison_cohort = comparison_group_ari,
    weight_formula = "covid~age_group+sex_cat+race_eth_cat+hospitalized_at_index+ce_week_factor+ highrisk_flag", ## PULMONARY_RESP?
    iptw_description = "resp_w_flu_site_exc_resp_15_highrisk")

dataset_iptw_ari <- results_tbl( paste0("resp_w_flu_site_exc_resp_15_highrisk",'_tabc')) %>%
  collect()%>%
  set_iptw_levels()
# mutate(hospitalized_at_index = factor(hospitalized_at_index, levels = c("Not hospitalized", "Hospitalized"))) %>%
# mutate(age_group = factor(age_group, levels = c("<6 months", "6m-1y old", "1-2y old", "2-3y old", "3-4y old", "4-5y old"))) %>%
# mutate(sex_cat = relelvel(factor(sex_cat), ref = "Female")) %>%


models_ari_2b <- run_logistic_regression(df = dataset_iptw_ari,
                                         model_formula = paste0(outcome_ari,
                                                                formula_no_pulmonary, " + highrisk_flag"),
                                         weight_method = "lr_s")

base_data_2b <- as.data.frame(coef(models_ari_2b))
base_data_2b <- cbind.data.frame(base_data_2b, confint(models_ari_2b))
base_data_2b <- cbind.data.frame(base_data_2b, as.data.frame(coef(summary((models_ari_2b)))[,4]))

var_to_omit <- c("(Intercept)", "util_inpatientno_visits", "util_inpatientlow_utilizer", "util_inpatienthigh_utilizer",
                 "util_inpatientmoderate_utilizer",
                 "util_outpatientno_visits", "util_outpatientlow_utilizer", "util_outpatienthigh_utilizer", "util_outpatientmoderate_utilizer",
                 "util_edno_visits", "util_edlow_utilizer", "util_edhigh_utilizer", "util_edmoderate_utilizer",
                 "util_otherno_visits", "util_otherlow_utilizer", "util_otherhigh_utilizer", "util_othermoderate_utilizer")

base_data_2b <- base_data_2b %>% rename(OR = "coef(models_ari_2b)", 
                                   pval = "coef(summary((models_ari_2b)))[, 4]", 
                                   lower = "2.5 %", upper = "97.5 %") %>%
  mutate(covariates = rownames(.)) %>%
  filter(!(covariates %in% var_to_omit)) %>%
  mutate(across(c("OR", "lower", "upper"), ~exp(.x))) %>%
  mutate(across(c("OR", "lower", "upper"), ~round(.x, 2))) %>%
  mutate(ci = paste0(lower, "-", upper)) %>%
  mutate(covariates = case_when(covariates == "highrisk_flagYes"~ "At least 1 high risk condition", 
                                covariates == "hospitalized_at_indexYes" ~ "Hopspitalized for index infection",
                                covariates == "race_eth_catOther/Unknown" ~ "Other/Unknown race",
                                covariates == "race_eth_catHispanic" ~ "Hispanic",
                                covariates == "race_eth_catNH_White" ~ "Non-Hispanic White",
                                covariates == "race_eth_catNH_Black/AA" ~ "Non-Hispanic Black/AA",
                                covariates == "race_eth_catNH_Asian/PI" ~ "Non-Hispanic Asian/PI",
                                covariates == "race_eth_catNH_Multiple" ~ "Non-Hispanic Multiple races",
                                covariates == "sex_catMale" ~ "Male",
                                covariates == "sex_catFemale" ~ "Female",
                                covariates == "age_group<6 months" ~ "<6 months",
                                covariates == "age_group6m-1y old" ~ "6m-1y old",
                                covariates == "age_group1-2y old" ~ "1-2y old",
                                covariates == "age_group2-3y old" ~ "2-3y old",
                                covariates == "age_group3-4y old" ~ "3-4y old",
                                covariates == "age_group4-5y old" ~ "4-5y old",
                                covariates == "covid" ~ "SARS-CoV-2 Infections",
                                TRUE ~ covariates))%>%
  mutate(pval = case_when(pval > 0.05 ~ format(round(pval, 4), scientific = FALSE),
                          pval >= 0.005 ~ paste0(format(round(pval, 4), scientific = FALSE), "*"),
                          pval >= 0.001 ~ paste0(format(round(pval, 4), scientific = FALSE), "**"),
                          pval >= 0.0001 ~ paste0(format(round(pval, 4), scientific = FALSE), "***"),
                          TRUE ~ "<0.0001****"))

png("results/figures_10082024/Figure_2b_ari15_180_highrisk_site_excl_w_adj.png", width = 20, height = 10, units = "cm", res = 300 )
p2b_w_adj <- base_data_2b|>
  forestplot(labeltext = c(covariates, OR, ci, pval),
             mean = OR, 
             clip = c(0.1, 5.0), #xlim
             xlog = FALSE,
             xlab = "Odds ratios", boxsize = 0.25, zero = 1) |>
  fp_set_style(box = "darkblue",
               line = "darkblue") |>
  fp_add_header(covariates = c("", ""),
                ci = c("95% CI", ""),
                OR = c("OR", ""),
                pval = c("p-value", ""))      
print(p2b_w_adj)

dev.off()

p_s3b <- ari_iptw_15_180$love_plot + 
  geom_vline(xintercept = 0.2, linetype = "dashed") +
  geom_vline(xintercept = -0.2, linetype = "dashed")
print(p_s3b)
ggsave("results/figures_10082024/Figure_s3b_ari_15_180_highrisk_site_excl_w_adj.png",p_s3b,  width = 6.5, height = 10, unit = "in")

## New Fig 2c: any infection 15-180, 1 model: weighted adjusted
outcome_any = "outcome_general15d"
comparison_group_any = "Influenza"

formula_no_pulmonary = "~ covid+ age_group + race_eth_cat + sex_cat + hospitalized_at_index +  util_inpatient +  util_outpatient +  util_ed + util_other"

any_iptw_15_180 <- analytic_dataset_final_resp %>%
  #mutate(ce_month = floor_date(ce_date, "month")) %>%
  #mutate(ce_month_s = as.character(ce_month)) %>%
  #filter(resp_outcome_is_reinfection=="0") %>% # TODO allow for re-infection
  mutate(resp15d_or_covid15d = case_when(outcome_resp15d==1 ~ 1,
                                         outcome_covid15d==1 ~ 1,
                                         TRUE ~ 0)) %>%
  # mutate(hospitalized_at_index = ifelse(hospitalized_at_index == "Yes", 1, 0)) %>%
  # mutate(highrisk_flag = ifelse(highrisk_flag == "Yes", 1, 0))%>%
  do_iptw(#analytic_dataset = ari_iptw_15_180, 
    cohort_entry_start_date = cohort_entry_start_date_resp,
    comparison_cohort = comparison_group_any,
    weight_formula = "covid~age_group+sex_cat+race_eth_cat+hospitalized_at_index+ce_week_factor+ highrisk_flag", ## PULMONARY_RESP?
    iptw_description = "resp_covid_flu_site_exc_general15d_highrisk")

dataset_iptw_any <- results_tbl( paste0("resp_covid_flu_site_exc_general15d_highrisk",'_tabc')) %>%
  collect()%>%
  set_iptw_levels()
# mutate(hospitalized_at_index = factor(hospitalized_at_index, levels = c("Not hospitalized", "Hospitalized"))) %>%
# mutate(age_group = factor(age_group, levels = c("<6 months", "6m-1y old", "1-2y old", "2-3y old", "3-4y old", "4-5y old"))) %>%
# mutate(sex_cat = relelvel(factor(sex_cat), ref = "Female")) %>%


models_any_2c <- run_logistic_regression(df = dataset_iptw_any,
                                         model_formula = paste0(outcome_any,
                                                                formula_no_pulmonary, " + highrisk_flag"),
                                         weight_method = "lr_s")

base_data_2c <- as.data.frame(coef(models_any_2c))
base_data_2c <- cbind.data.frame(base_data_2c, confint(models_any_2c))
base_data_2c <- cbind.data.frame(base_data_2c, as.data.frame(coef(summary((models_any_2c)))[,4]))

var_to_omit <- c("(Intercept)", "util_inpatientno_visits", "util_inpatientlow_utilizer", "util_inpatienthigh_utilizer",
                 "util_inpatientmoderate_utilizer",
                 "util_outpatientno_visits", "util_outpatientlow_utilizer", "util_outpatienthigh_utilizer", "util_outpatientmoderate_utilizer",
                 "util_edno_visits", "util_edlow_utilizer", "util_edhigh_utilizer", "util_edmoderate_utilizer",
                 "util_otherno_visits", "util_otherlow_utilizer", "util_otherhigh_utilizer", "util_othermoderate_utilizer")

base_data_2c <- base_data_2c %>% rename(OR = "coef(models_any_2c)", 
                                        pval = "coef(summary((models_any_2c)))[, 4]", 
                                        lower = "2.5 %", upper = "97.5 %") %>%
  mutate(covariates = rownames(.)) %>%
  filter(!(covariates %in% var_to_omit)) %>%
  mutate(across(c("OR", "lower", "upper"), ~exp(.x))) %>%
  mutate(across(c("OR", "lower", "upper"), ~round(.x, 2))) %>%
  mutate(ci = paste0(lower, "-", upper)) %>%
  mutate(covariates = case_when(covariates == "highrisk_flagYes"~ "At least 1 high risk condition", 
                                covariates == "hospitalized_at_indexYes" ~ "Hopspitalized for index infection",
                                covariates == "race_eth_catOther/Unknown" ~ "Other/Unknown race",
                                covariates == "race_eth_catHispanic" ~ "Hispanic",
                                covariates == "race_eth_catNH_White" ~ "Non-Hispanic White",
                                covariates == "race_eth_catNH_Black/AA" ~ "Non-Hispanic Black/AA",
                                covariates == "race_eth_catNH_Asian/PI" ~ "Non-Hispanic Asian/PI",
                                covariates == "race_eth_catNH_Multiple" ~ "Non-Hispanic Multiple races",
                                covariates == "sex_catMale" ~ "Male",
                                covariates == "sex_catFemale" ~ "Female",
                                covariates == "age_group<6 months" ~ "<6 months",
                                covariates == "age_group6m-1y old" ~ "6m-1y old",
                                covariates == "age_group1-2y old" ~ "1-2y old",
                                covariates == "age_group2-3y old" ~ "2-3y old",
                                covariates == "age_group3-4y old" ~ "3-4y old",
                                covariates == "age_group4-5y old" ~ "4-5y old",
                                covariates == "covid" ~ "SARS-CoV-2 Infections",
                                TRUE ~ covariates))%>%
  mutate(pval = case_when(pval > 0.05 ~ format(round(pval, 4), scientific = FALSE),
                          pval >= 0.005 ~ paste0(format(round(pval, 4), scientific = FALSE), "*"),
                          pval >= 0.001 ~ paste0(format(round(pval, 4), scientific = FALSE), "**"),
                          pval >= 0.0001 ~ paste0(format(round(pval, 4), scientific = FALSE), "***"),
                          TRUE ~ "<0.0001****"))

png("results/figures_10082024/Figure_2c_any15_180_highrisk_site_excl_w_adj.png", width = 20, height = 10, units = "cm", res = 300 )
p2c_w_adj <- base_data_2c|>
  forestplot(labeltext = c(covariates, OR, ci, pval),
             mean = OR, 
             clip = c(0.1, 5.0), #xlim
             xlog = FALSE,
             xlab = "Odds ratios", boxsize = 0.25, zero = 1) |>
  fp_set_style(box = "darkblue",
               line = "darkblue") |>
  fp_add_header(covariates = c("", ""),
                ci = c("95% CI", ""),
                OR = c("OR", ""),
                pval = c("p-value", ""))      
print(p2c_w_adj)

dev.off()

p_s3c <- any_iptw_15_180$love_plot + 
  geom_vline(xintercept = 0.2, linetype = "dashed") +
  geom_vline(xintercept = -0.2, linetype = "dashed")
print(p_s3c)
ggsave("results/figures_10082024/Figure_s3c_any_15_180_highrisk_site_excl_w_adj.png",p_s3c,  width = 6.5, height = 10, unit = "in")


############################## Suchitra also wants a figure 2a to compare covid vs respiratory cohort ###################

outcome_15_180_resp = "outcome_rsv_15_to_180"
comparison_group_resp = "Respiratory"

formula_no_pulmonary = "~ covid+ age_group + race_eth_cat + sex_cat + hospitalized_at_index +  util_inpatient +  util_outpatient +  util_ed + util_other"

rsv_iptw_15_180_highrisk <- analytic_dataset_final %>%
  do_iptw(cohort_entry_start_date = cohort_entry_start_date,
          comparison_cohort = comparison_group_resp,
          weight_formula = "covid~age_group+sex_cat+race_eth_cat+ce_week_factor",
          iptw_description = "covid_resp_rsv15_180_highrisk_site_excl")

df_15_180_highrisk <- results_tbl( paste0("covid_resp_rsv15_180_highrisk_site_excl",'_tabc')) %>%
  collect() %>% 
  mutate(hospitalized_at_index = factor(hospitalized_at_index, levels = c("Not hospitalized", "Hospitalized"))) %>%
  mutate(age_group = factor(age_group, levels = c("<6 months", "6m-1y old", "1-2y old", "2-3y old", "3-4y old", "4-5y old"))) %>%
  mutate(sex_cat = relevel(factor(sex_cat), ref = "Male"))

# weighted, adjusted models
models_highrisk_5 <- run_logistic_regression(df = df_15_180_highrisk,
                                              model_formula = paste0(outcome_15_180, formula_no_pulmonary, "+ highrisk_flag"),
                                              weight_method = "lr")

base_data_5 <- as.data.frame(coef(models_highrisk_5))
base_data_5 <- cbind.data.frame(base_data_5, confint(models_highrisk_5))
base_data_5 <- cbind.data.frame(base_data_5, as.data.frame(coef(summary(models_highrisk_5))[,4]))

var_to_omit <- c("(Intercept)", "util_inpatientno_visits", "util_inpatientlow_utilizer", "util_inpatienthigh_utilizer",
                 "util_inpatientmoderate_utilizer",
                 "util_outpatientno_visits", "util_outpatientlow_utilizer", "util_outpatienthigh_utilizer", "util_outpatientmoderate_utilizer",
                 "util_edno_visits", "util_edlow_utilizer", "util_edhigh_utilizer", "util_edmoderate_utilizer",
                 "util_otherno_visits", "util_otherlow_utilizer", "util_otherhigh_utilizer", "util_othermoderate_utilizer")

base_data_5 <- base_data_5 %>% rename(OR = "coef(models_highrisk_5)", 
                                        pval = "coef(summary(models_highrisk_5))[, 4]", 
                                        lower = "2.5 %", upper = "97.5 %") %>%
  mutate(covariates = rownames(.)) %>%
  filter(!(covariates %in% var_to_omit)) %>%
  mutate(across(c("OR", "lower", "upper"), ~exp(.x))) %>%
  mutate(across(c("OR", "lower", "upper"), ~round(.x, 2))) %>%
  mutate(ci = paste0(lower, "-", upper)) %>%
  mutate(covariates = case_when(covariates == "highrisk_flagYes"~ "At least 1 high risk condition", 
                                covariates == "hospitalized_at_indexHospitalized" ~ "Hopspitalized for index infection",
                                covariates == "race_eth_catOther/Unknown" ~ "Other/Unknown race",
                                covariates == "race_eth_catHispanic" ~ "Hispanic",
                                covariates == "race_eth_catNH_White" ~ "Non-Hispanic White",
                                covariates == "race_eth_catNH_Black/AA" ~ "Non-Hispanic Black/AA",
                                covariates == "race_eth_catNH_Asian/PI" ~ "Non-Hispanic Asian/PI",
                                covariates == "race_eth_catNH_Multiple" ~ "Non-Hispanic Multiple races",
                                covariates == "sex_catMale" ~ "Male",
                                covariates == "sex_catFemale" ~ "Female",
                                covariates == "age_group<6 months" ~ "<6 months",
                                covariates == "age_group6m-1y old" ~ "6m-1y old",
                                covariates == "age_group1-2y old" ~ "1-2y old",
                                covariates == "age_group2-3y old" ~ "2-3y old",
                                covariates == "age_group3-4y old" ~ "3-4y old",
                                covariates == "age_group4-5y old" ~ "4-5y old",
                                covariates == "covid" ~ "SARS-CoV-2 Infections",
                                TRUE ~ covariates))%>%
  mutate(pval = case_when(pval > 0.05 ~ format(round(pval, 4), scientific = FALSE),
                          pval >= 0.005 ~ paste0(format(round(pval, 4), scientific = FALSE), "*"),
                          pval >= 0.001 ~ paste0(format(round(pval, 4), scientific = FALSE), "**"),
                          pval >= 0.0001 ~ paste0(format(round(pval, 4), scientific = FALSE), "***"),
                          TRUE ~ "<0.0001****"))

png("results/figures_10082024/Figure_5_rsv15_180_highrisk_site_excl_w_adj.png", width = 20, height = 10, units = "cm", res = 300 )
p5_w_adj <- base_data_5 |>
  forestplot(labeltext = c(covariates, OR, ci, pval),
             mean = OR, 
             clip = c(0.1, 5.0), #xlim
             xlog = FALSE,
             xlab = "Odds ratios", boxsize = 0.25, zero = 1) |>
  fp_set_style(box = "darkblue",
               line = "darkblue") |>
  fp_add_header(covariates = c("", ""),
                ci = c("95% CI", ""),
                OR = c("OR", ""),
                pval = c("p-value", ""))      
print(p5_w_adj)

dev.off()
