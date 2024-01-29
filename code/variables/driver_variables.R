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
results_tbl("base_cohort") %>% left_join(cdm_tbl("visit_occurrence") %>% select(visit_occurrence_id, person_id, visit_start_date), by=c("person_id", "ce_date"="visit_start_date")) %>% left_join(cdm_tbl("visit_payer") %>% select(visit_occurrence_id, plan_class), by="visit_occurrence_id") %>% group_by(plan_class) %>% summarise(n())

