### File for actual creation of the cohort

## Variables that will be needed:
## covid evidence or influenza evidence
## utilization visits before/after
## washouts applied: serology and covid test
## exclude co-infection or flag
## age applied
### TODO remove file, obsolete now 

### Cohort data structure:

## cols:
# columns <- c("person_id", 
#              "ce_date",
#              "sub_cohort" # Formerly cohort_assignment
#              )
# 
# ## Minimum number of columns needed to create the base cohort, and its tidy
# ## TODO: need to maintain visit_occurrence_id for the index visit for each person, otherwise won't know which to use for visit_payer crosswalk? could reverse engineer
# 
# ### Washouts
# ce_start_date = as.Date("2021-04-01")
# ce_end_date = as.Date("2023-10-01")
# # cohort_1_label = "_06210223_"
# cohort_1_label = "_04211023_"
# odr <- cdm_tbl("observation_derivation_recover") %>% 
#   filter(observation_date < ce_end_date, observation_date >= ce_start_date)
# 
# 
# ## Code continues from `driver_cohort_attrition`
# 
# cohort_with_washouts <-
#   results_tbl(paste0("ce", cohort_1_label, "unfiltered_demo")) %>% 
#   mutate(sub_cohort = cohort_assignment) %>% 
#   apply_washout_logic(odr = odr) %>% 
#   compute_new()
# 
# # cohort_with_washouts %>% 
# #   group_by(study_eligible, washout) %>% 
# #   summarise(n())
# 
# # cohort_with_washouts %>% 
# #   filter(cohort_assignment %in% c("Covid", "Influenza")) %>% 
# #   summarise_as_fraction(group = c(washout, study_eligible), distinct_key = person_id)
# 
# cohort_with_washouts %>% 
#   filter(study_eligible == 1) %>% 
#   filter(washout == 0) %>% 
#   filter(cohort_assignment != "No infection") %>% 
#   select(person_id, ce_date, sub_cohort = cohort_assignment) %>% 
#   output_tbl("base_cohort", indexes=c("person_id"))
# 
# 
# 
# 
# 
# 
