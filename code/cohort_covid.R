flag_covid_positives <- function(cohort, odr_tbl) {
  cohort %>% 
    flag_positive_sars_cov2_test(odr_tbl) %>% 
    flag_code_u099(odr_tbl) %>%
    # flag_prescription_pax_des(odr_tbl) %>% # add scrip_rx_date if applying
    flag_code_u07_one_or_more(odr_tbl) %>%
    # flag_serology(odr_tbl) %>% 
    mutate(earliest_cov_event_date = pmin(test_date, covid_dx_date, pasc_dx_date)) %>% 
    impute_pasc_orphans() %>% 
    mutate(covid_event_positive = ifelse(!is.na(earliest_cov_event_date), 1, 0)) %>% 
    return()
}

flag_positive_sars_cov2_test <- function(cohort, odr_tbl) {
  # Logic to add a flag if a patient has a positive sars cov 2 test
  # Just grabbing earliest covid test date, per person_id
  
  pos_test_cohort <- odr_tbl %>% 
    filter(observation_concept_id %in% c(2000001530, 2000001529), value_as_concept_id %in% c(9191L, 2000001526)) %>% 
    select(person_id, test_date = observation_date) %>% 
    group_by(person_id) %>%
    slice_min(test_date, with_ties = FALSE) %>%
    mutate(positive_test_flag = 1) %>% 
    compute_new(indexes=c("person_id"))
  
  cohort %>% 
    left_join(
      pos_test_cohort,
      by=c("person_id")
    ) %>% 
    mutate(positive_test_flag = ifelse(is.na(positive_test_flag), 0, positive_test_flag)) %>% 
    return()
}

flag_code_u07_one_or_more <- function(cohort, odr_tbl) {
  # Logic to add a flag for a cohort having one or mor U07.1 codes (covid diagnosis)
  covid_dx_tbl <- odr_tbl %>% 
    filter(observation_concept_id == 2000001527, value_as_concept_id == 2000001525) %>% 
    select(person_id, covid_dx_date = observation_date) %>%
    group_by(person_id) %>%
    slice_min(covid_dx_date, with_ties = FALSE) %>%
    mutate(covid_dx_flag = 1) %>% 
    compute_new(indexes=c("person_id"))
  
  cohort %>% 
    left_join(covid_dx_tbl, by="person_id") %>% 
    mutate(covid_dx_flag = ifelse(is.na(covid_dx_flag), 0, covid_dx_flag)) %>% 
    return()
}

flag_code_u099 <- function(cohort, odr_tbl) {
  # Logic to add a flag if a patient has a u09 dx code
  pasc_flag <- odr_tbl %>%
    filter(observation_concept_id == 2000001527, value_as_concept_id %in% c(2000001520, 2000001533, 703578)) %>%
    select(person_id, pasc_dx_date = observation_date) %>%
    group_by(person_id) %>%
    slice_min(pasc_dx_date, with_ties = FALSE) %>%
    mutate(pasc_dx_flag = 1) %>% 
    compute_new(indexes=c("person_id"))
  
  cohort %>%
    left_join(pasc_flag, by="person_id") %>% 
    mutate(pasc_dx_flag = ifelse(is.na(pasc_dx_flag), 0, pasc_dx_flag)) %>% 
    return()
}

## Takes the cohort with the following columns, at least:
## test_date, covid_dx_date, pasc_dx_date, scrip_rx_date
## positive_test_flag, covid_dx_flag, scrip_rx_flag, pasc_dx_flag
## IF the other dates and flags are null/0, and ONLY the pasc_dx_flag is true,
## that means the PASC is the earliest covid event, therefore, the ce_date needs to be imputed properly.
## Will just add a new column that is index date (same as ce_date) except for the u09 patients
impute_pasc_orphans <- function(cohort) {
  ## 
  cohort %>% 
    mutate(index_date = case_when(
      positive_test_flag == 1 |
        covid_dx_flag == 1 ~ earliest_cov_event_date,
      pasc_dx_flag == 1 &
        (positive_test_flag == 0 &
           covid_dx_flag == 0 ) ~ as.Date(pasc_dx_date - days(59))
    )) %>% 
    mutate(index_date_imputed = case_when(
      positive_test_flag == 1 |
        covid_dx_flag == 1 ~ 0,
      pasc_dx_flag == 1 &
        positive_test_flag == 0 &
        covid_dx_flag == 0 ~ 1)) %>% 
    return()
}