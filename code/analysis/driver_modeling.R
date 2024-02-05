### Driver file for taking up an analytic dataset with outcomes and covariates, performing IPW weighting and Cox modeling
## for results

library(survey)
library(survival)
library(WeightIt)

### Making it into an analytic Cox format dataset
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

create_cox_format_df <- function(outcome_data) {
  cox_df <- outcome_data %>% 
    filter(sub_cohort %in% c("Covid", "Influenza")) %>% 
    mutate(exposure = case_when(sub_cohort=="Covid" ~ 1,
                                sub_cohort == "Influenza" ~ 0)) %>% 
    mutate(ce_month = as.Date(floor_date(ce_date, unit="month"))) %>% 
    mutate(ce_days_secular = as.numeric(ce_month - as.Date("2021-06-01"))) %>% 
    mutate(days_til_earliest_resp_outcome = ifelse(is.na(days_til_earliest_resp_outcome), 180, days_til_earliest_resp_outcome)) %>% 
    mutate(days_til_earliest_general_outcome = ifelse(is.na(days_til_earliest_general_outcome), 180, days_til_earliest_general_outcome)) %>% 
    mutate(days_til_earliest_rsv_outcome = ifelse(is.na(days_til_earliest_rsv_outcome), 180, days_til_earliest_rsv_outcome)) %>% 
    compute_new()
  
  cox_df %>% 
    return()
}

add_iptw_weights <- function(cox_format_df, estimand, formula_expression, output_tbl_name) {
  
  cox_df_collected <- cox_format_df %>% collect()
  weights.out <- weightit(as.formula(formula_expression), data=cox_df_collected,
                          method="glm", estimand=estimand) #focal = "circle pill")
  cox_df_collected$wts = weights.out$weights
  
  ## TODO
  # Need to figure out the best way to output table metadata on the IPTW balance plots and SMDs
  weight_summary <- summary(weights.out)
  
  weights_metadata = list()
  weights_metadata$summary <- weight_summary
  
  cox_df_collected %>% 
    output_tbl(output_tbl_name)
  
  weights_metadata %>% 
    return()
  
}

weights.out <- weightit(as.formula("exposure~ce_days_secular"), data=cohort_outcomes,
                        method="glm", estimand="ATE") #focal = "circle pill")
cohort_outcomes$w = weights.out$weights

summary(weights.out)

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

perform_cox_model <- function(cox_data, outcome, time_until) {
  cox_mod_adj <- coxph(Surv(days_til_earliest_general_outcome, has_postacute_general_outcome) ~ 
                         exposure + ce_days_secular + race_eth_cat, #need the time to event TODO start here on Tuesday
                       data=cox_data,
                       weights = wts)
  
  summary(cox_mod_adj) %>% 
    return()
}

SurvFun <- function(fun.time, fun.event, grouping = 1, fun.dat) {
  params <- list(fun.time = substitute(fun.time),
                 fun.event = substitute(fun.event),
                 grouping = substitute(grouping), 
                 fun.dat = substitute(fun.dat))
  expr <- substitute(survfit(Surv(time = fun.time, event = fun.event) ~ grouping, 
                             data = fun.dat), params)
  eval.parent(expr)
}




