### File for actual creation of the cohort

## Variables that will be needed:
## covid evidence or influenza evidence
## utilization visits before/after
## washouts applied: serology and covid test
## exclude co-infection or flag
## age applied

### Cohort data structure:

## cols:
columns <- c("person_id", 
             "ce_date",
             "sub_cohort" # Formerly cohort_assignment
             )

## Minimum number of columns needed to create the base cohort, and its tidy
## TODO: need to maintain visit_occurrence_id for the index visit for each person, otherwise won't know which to use for visit_payer crosswalk? could reverse engineer

### Washouts
ce_start_date = as.Date("2021-04-01")
ce_end_date = as.Date("2023-10-01")
# cohort_1_label = "_06210223_"
cohort_1_label = "_04211023_"
odr <- cdm_tbl("observation_derivation_recover") %>% 
  filter(observation_date < ce_end_date, observation_date >= ce_start_date)

### Performing the washouts
## Getting data pertaining to the washouts
#MISC or PASC: if observation_concept_id==2000001527, value_as_concept_id==703578, or 2000001520  and patient does NOT have covid, then washout
# Serology test, observation_concept_id==2000001528, value_as_concept_id %in% c(2000001526, 9191)
apply_washout_logic <- function(cohort, odr) {
  cohort %>% 
    left_join(odr %>% 
                filter(observation_concept_id %in% c(2000001527, 2000001528),
                       value_as_concept_id %in% c(703578, 2000001520, 2000001526, 9191)) %>% 
                select(person_id, observation_date, observation_concept_id, value_as_concept_id),
              by="person_id") %>% 
    mutate(washout_reason = case_when(
      observation_concept_id==2000001527 & 
        value_as_concept_id %in% c(703578, 2000001520) & 
        sub_cohort %in% c("Influenza", "No infection") ~ "non_covid pt with_misc or pasc",
      observation_concept_id==2000001528 & 
        value_as_concept_id %in% c(2000001526, 9191) & 
        sub_cohort == "Covid" &
        observation_date < ce_date ~ "covid pat with serology pos prior to ce date",
      observation_concept_id==2000001528 & 
        value_as_concept_id %in% c(2000001526, 9191) & 
        sub_cohort %in% c("Influenza", "No infection") ~ "non_covid pat with serology pos",
      TRUE ~ "none"
    )) %>% 
    mutate(washout = case_when(
      washout_reason == "none" ~ 0,
      TRUE ~ 1
    )) %>% 
    group_by(person_id) %>% 
    slice_max(washout, with_ties = FALSE) %>% 
    ungroup() %>% 
    return()
}


## Code continues from `driver_cohort_attrition`

cohort_with_washouts <-
  results_tbl(paste0("ce", cohort_1_label, "unfiltered")) %>% 
  mutate(sub_cohort = cohort_assignment) %>% 
  apply_washout_logic(odr = odr) %>% 
  compute_new()

cohort_with_washouts %>% 
  group_by(study_eligible, washout) %>% 
  summarise(n())

cohort_with_washouts %>% 
  filter(cohort_assignment %in% c("Covid", "Influenza")) %>% 
  summarise_as_fraction(group = c(washout, study_eligible), distinct_key = person_id)

cohort_with_washouts %>% 
  filter(study_eligible == 1) %>% 
  filter(washout == 0) %>% 
  filter(cohort_assignment != "No infection") %>% 
  select(person_id, ce_date, sub_cohort = cohort_assignment) %>% 
  output_tbl("base_cohort", indexes=c("person_id"))






