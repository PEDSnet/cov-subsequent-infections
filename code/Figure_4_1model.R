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

formula_no_pulmonary = "~ covid+ age_group + race_eth_cat + sex_cat + hospitalized_at_index +  util_inpatient +  util_outpatient +  util_ed + util_other"

rsv_iptw_15_300_highrisk <- analytic_dataset_final %>%
  do_iptw(cohort_entry_start_date = cohort_entry_start_date,
          comparison_cohort = comparison_group,
          weight_formula = "covid~age_group+sex_cat+race_eth_cat+ce_week_factor",
          iptw_description = "covid_flu_rsv15_300_highrisk_site_excl")

df_15_300_highrisk <- results_tbl( paste0("covid_flu_rsv15_300_highrisk_site_excl",'_tabc')) %>%
  collect() %>% 
  mutate(hospitalized_at_index = factor(hospitalized_at_index, levels = c("Not hospitalized", "Hospitalized"))) %>%
  mutate(age_group = factor(age_group, levels = c("<6 months", "6m-1y old", "1-2y old", "2-3y old", "3-4y old", "4-5y old"))) %>%
  mutate(sex_cat = relevel(factor(sex_cat), ref = "Male"))

# weighted, adjusted models
models_highrisk_4 <- run_logistic_regression(df = df_15_300_highrisk,
                                              model_formula = paste0(outcome_15_300, formula_no_pulmonary, "+ highrisk_flag"),
                                              weight_method = "lr")

base_data_4 <- as.data.frame(coef(models_highrisk_4))
base_data_4 <- cbind.data.frame(base_data_4, confint(models_highrisk_4))
base_data_4 <- cbind.data.frame(base_data_4, as.data.frame(coef(summary(models_highrisk_4))[,4]))

var_to_omit <- c("(Intercept)", "util_inpatientno_visits", "util_inpatientlow_utilizer", "util_inpatienthigh_utilizer",
                 "util_inpatientmoderate_utilizer",
                 "util_outpatientno_visits", "util_outpatientlow_utilizer", "util_outpatienthigh_utilizer", "util_outpatientmoderate_utilizer",
                 "util_edno_visits", "util_edlow_utilizer", "util_edhigh_utilizer", "util_edmoderate_utilizer",
                 "util_otherno_visits", "util_otherlow_utilizer", "util_otherhigh_utilizer", "util_othermoderate_utilizer")

base_data_4a <- base_data_4 %>% rename(OR = "coef(models_highrisk_4)", 
                                        pval = "coef(summary(models_highrisk_4))[, 4]", 
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

png("results/figures_10082024/Figure_4a_rsv15_300_highrisk_site_excl_w_adj.png", width = 20, height = 10, units = "cm", res = 300 )
p4a_w_adj <- base_data_4a |>
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
print(p4a_w_adj)

dev.off()


message("Figure 4b: RSV within 30-180 days (past post-acute)")
outcome_30_180 = "outcome_rsv_30_to_180"
comparison_group = "Influenza"
rsv_iptw_30_180_highrisk <- analytic_dataset_final %>%
  do_iptw(cohort_entry_start_date = cohort_entry_start_date,
          comparison_cohort = comparison_group,
          weight_formula = "covid~age_group+sex_cat+race_eth_cat+ce_week_factor",
          iptw_description = "covid_flu_rsv30_180_highrisk_site_excl")

df_30_180_highrisk <- results_tbl( paste0("covid_flu_rsv30_180_highrisk_site_excl",'_tabc')) %>%
  collect() %>% 
  mutate(hospitalized_at_index = factor(hospitalized_at_index, levels = c("Not hospitalized", "Hospitalized"))) %>%
  mutate(age_group = factor(age_group, levels = c("<6 months", "6m-1y old", "1-2y old", "2-3y old", "3-4y old", "4-5y old"))) %>%
  mutate(sex_cat = relevel(factor(sex_cat), ref = "Male"))

# weighted, adjusted models
models_highrisk_4a <- run_logistic_regression(df = df_30_180_highrisk,
                                             model_formula = paste0(outcome_30_180, formula_no_pulmonary, "+ highrisk_flag"),
                                             weight_method = "lr")

base_data_4a <- as.data.frame(coef(models_highrisk_4a))
base_data_4a <- cbind.data.frame(base_data_4a, confint(models_highrisk_4a))
base_data_4a <- cbind.data.frame(base_data_4a, as.data.frame(coef(summary(models_highrisk_4a))[,4]))

var_to_omit <- c("(Intercept)", "util_inpatientno_visits", "util_inpatientlow_utilizer", "util_inpatienthigh_utilizer",
                 "util_inpatientmoderate_utilizer",
                 "util_outpatientno_visits", "util_outpatientlow_utilizer", "util_outpatienthigh_utilizer", "util_outpatientmoderate_utilizer",
                 "util_edno_visits", "util_edlow_utilizer", "util_edhigh_utilizer", "util_edmoderate_utilizer",
                 "util_otherno_visits", "util_otherlow_utilizer", "util_otherhigh_utilizer", "util_othermoderate_utilizer")

base_data_4a <- base_data_4a %>% rename(OR = "coef(models_highrisk_4a)", 
                                       pval = "coef(summary(models_highrisk_4a))[, 4]", 
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

png("results/figures_10082024/Figure_4a_rsv_30_180_highrisk_site_excl_w_adj.png", width = 20, height = 10, units = "cm", res = 300 )
p4a_w_adj <- base_data_4a |>
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
print(p4a_w_adj)

dev.off()

message("Figure 4a: RSV within 0-180 days (include co-infections)")
outcome_0_180 = "outcome_rsv_0_to_180"
comparison_group = "Influenza"
# create this outcome from 0_to_60 and 61_to_180 outcomes
analytic_dataset_final <- analytic_dataset_final%>%
  mutate(outcome_rsv_0_to_180 = ifelse(outcome_rsv_0_to_60 == 1 | outcome_rsv_61_to_180 == 1, 1, 0)) %>%
  do_iptw(cohort_entry_start_date = cohort_entry_start_date,
          comparison_cohort = comparison_group,
          weight_formula = "covid~age_group+sex_cat+race_eth_cat+ce_week_factor+highrisk_flag",
          iptw_description = "covid_flu_outcome_rsv0_180_highrisk_site_excl")

df_0_180_highrisk <- results_tbl( paste0("covid_flu_outcome_rsv0_180_highrisk_site_excl",'_tabc')) %>%
  collect() %>% 
  mutate(hospitalized_at_index = factor(hospitalized_at_index, levels = c("Not hospitalized", "Hospitalized"))) %>%
  mutate(age_group = factor(age_group, levels = c("<6 months", "6m-1y old", "1-2y old", "2-3y old", "3-4y old", "4-5y old"))) %>%
  mutate(sex_cat = relevel(factor(sex_cat), ref = "Male"))

# weighted, adjusted models
models_highrisk_4b <- run_logistic_regression(df = df_0_180_highrisk,
                                              model_formula = paste0(outcome_0_180, formula_no_pulmonary, "+ highrisk_flag"),
                                              weight_method = "lr")

base_data_4b <- as.data.frame(coef(models_highrisk_4b))
base_data_4b <- cbind.data.frame(base_data_4b, confint(models_highrisk_4b))
base_data_4b <- cbind.data.frame(base_data_4b, as.data.frame(coef(summary(models_highrisk_4b))[,4]))

var_to_omit <- c("(Intercept)", "util_inpatientno_visits", "util_inpatientlow_utilizer", "util_inpatienthigh_utilizer",
                 "util_inpatientmoderate_utilizer",
                 "util_outpatientno_visits", "util_outpatientlow_utilizer", "util_outpatienthigh_utilizer", "util_outpatientmoderate_utilizer",
                 "util_edno_visits", "util_edlow_utilizer", "util_edhigh_utilizer", "util_edmoderate_utilizer",
                 "util_otherno_visits", "util_otherlow_utilizer", "util_otherhigh_utilizer", "util_othermoderate_utilizer")

base_data_4b <- base_data_4b %>% rename(OR = "coef(models_highrisk_4b)", 
                                        pval = "coef(summary(models_highrisk_4b))[, 4]", 
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

png("results/figures_10082024/Figure_4b_rsv_0_180_highrisk_site_excl_w_adj.png", width = 20, height = 10, units = "cm", res = 300 )
p4b_w_adj <- base_data_4b |>
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
print(p4b_w_adj)

dev.off()
