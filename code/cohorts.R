#'
#' This file contains functions to identify cohorts in a study.  Its
#' contents, as well as the contents of any R files whose name begins
#' with "cohort_", will be automatically sourced when the request is
#' executed.
#'
#' For simpler requests, it will often be possible to do all cohort
#' manipulation in functions within this file.  For more complex
#' requests, it may be more readable to separate different cohorts or
#' operations on cohorts into logical groups in different files.
#'

get_patients_under_age_limit <- function(person_tbl, visit_tbl, ce_start_date, ce_end_date, age_limit_years_rough) {
  # Gets all patients under 6 during the whole study entry period
  visit_tbl %>% 
    filter(visit_start_date >= ce_start_date, visit_start_date < ce_end_date) %>% 
    select(person_id, visit_start_date) %>% 
    group_by(person_id) %>% 
    slice_sample(n=1) %>% # random visit in the period ? might pose age restriction IF patient gets filtered out for being over 5 but they have an index date before, so make it plus one year
    ungroup() %>% 
    compute_new(indexes=c("person_id")) %>% 
    inner_join(person_tbl %>% 
                 select(person_id, birth_date), by="person_id") %>% 
    mutate(random_visit_age = as.numeric(visit_start_date-birth_date)) %>% 
    filter(random_visit_age < age_limit_years_rough*365.25, random_visit_age >= 0) %>% 
    compute_new(indexes=c("person_id")) %>% 
    return()
}

get_influenza_evidence <- function(cohort_tbl, ce_start_date, ce_end_date) {
  # flu_dates <- cohort_tbl %>% 
  #   mutate(flu_date = case_when(person_id %% 5 == 0 ~ as.Date(visit_start_date + days(11)),
  #                               person_id %% 3 == 0 ~ as.Date(visit_start_date + days(17)),
  #                               TRUE ~ NA))# arbitrary random index date in place of
  # 
  # flu_dates %>% 
  #   return()
  ## Return the cohort with a flu evidence within the study period
  ## columns: person_id, flu_date
  
  influenza_positives <-
    cohort_tbl %>%
    get_flu_evidence(min_date = ce_start_date, max_date = ce_end_date)
  
  cohort_tbl %>% 
    left_join(influenza_positives %>% 
                select(person_id, 
                       flu_date = earliest_flu_evidence,
                       flu_dx_date,
                       earliest_flu_test
                       ),
              by="person_id") %>% 
    return()
  ## TODO finish return type here 
}

get_covid_evidence <- function(cohort_tbl, odr_tbl) {
  # covid_dates <- cohort_tbl %>% 
  #   mutate(covid_date = case_when(person_id %% 7 == 0 ~ as.Date(visit_start_date + days(11)),
  #                               person_id %% 4 == 0 ~ as.Date(visit_start_date + days(3)),
  #                               person_id %% 13 == 0 ~ as.Date(visit_start_date + days(18)),
  #                               TRUE ~ NA))# arbitrary random index date in place of
  # 
  ## Return the cohort with a covid evidence within the study period
  ## columns: person_id, covid_date
  covid_positives <- 
    cohort_tbl %>% 
    flag_covid_positives(odr_tbl)
  
  cohort_tbl %>% 
    left_join(covid_positives %>% 
                select(person_id, covid_date = index_date,
                       covid_index_date_imputed = index_date_imputed,
                       cov_test_date = test_date, covid_dx_date, pasc_dx_date),
              by="person_id") %>% 
    return()
  
}

### This will be the primary use function
identify_cohort_group <- function(cohort_tbl, odr_tbl, ce_start_date, ce_end_date) {
  # Logic in here to decide what to do with kids who had both within the study period, or who had one or the other or none
  relevant_infections <-
    cohort_tbl %>% 
    get_covid_evidence(odr_tbl) %>% 
    get_influenza_evidence(ce_start_date, ce_end_date) %>% 
    compute_new(indexes=c("person_id"))
  
  ## Apply rules:
  relevant_infections %>% 
    mutate(cohort_assignment = case_when(
      is.na(covid_date) & is.na(flu_date) ~ "No infection",
      is.na(covid_date) & !is.na(flu_date) ~ "Influenza",
      is.na(flu_date) & !is.na(covid_date) ~ "Covid",
      !is.na(flu_date) & !is.na(covid_date) & covid_date < flu_date ~ "Covid prior to influenza",
      !is.na(flu_date) & !is.na(covid_date) & flu_date < covid_date ~ "Influenza prior to covid",
      !is.na(flu_date) & !is.na(covid_date) & covid_date == flu_date ~ "Covid and influenza same day",
      TRUE ~ "Undiscernable"
    )) %>% 
    mutate(ce_date = case_when(
      cohort_assignment == "No infection" ~ visit_start_date, # later make random, could also pick visit as random from the start
      cohort_assignment == "Influenza" ~ flu_date,
      cohort_assignment == "Covid" ~ covid_date,
      cohort_assignment == "Covid prior to influenza" ~ covid_date,
      cohort_assignment == "Influenza prior to covid" ~ flu_date,
      cohort_assignment == "Covid and influenza same day" ~ covid_date,
      TRUE ~ NA
    )) %>% 
    select(-visit_start_date) %>% 
    mutate(age_years_on_ce_date = as.numeric(ce_date-birth_date)/365.25) %>% 
    return()
  
  
}

flag_study_eligiblity <- function(cohort_with_flags, ce_start_date, ce_end_date) {
  cohort_with_flags %>% 
    mutate(study_eligible = case_when(
      age_at_ce_under_limit==1 &
        visit_within_n_years>0 &
        has_visit_in_followup_period==1 &
        ce_date >= ce_start_date &
        ce_date < ce_end_date ~ 1,
      TRUE ~ 0
    )) %>% 
    return()
}

flag_utilization_level <- function(cohort, visit_tbl, lookback_days) {
  cohort_visits_lookback <- cohort %>% 
    select(person_id, ce_date) %>% 
    left_join(visit_tbl, by="person_id") %>% 
    filter(visit_start_date >= as.Date(ce_date - days(lookback_days)),
           visit_start_date < ce_date) %>% 
    compute_new(indexes=c("person_id"))
  
  ### What do we want? A column for each visit type, with the quintile category that that patient belongs to
  
  visit_lookback_counts <-
    cohort_visits_lookback %>% 
    mutate(visit_type_rough = case_when(visit_concept_name %in% c("Outpatient Visit - Non Physician",
                                                                  "Outpatient Visit") ~ "outpatient",
                                        visit_concept_name %in% c("Inpatient Visit - Ongoing", 
                                                                  "Inpatient Visit",
                                                                  "Observation Stay - PCORNet",
                                                                  "Emergency Department Admit to Inpatient Hospital Stay- PCORNet",
                                                                  "Emergency Room Visit") ~ "inpatient_or_ed",
                                        visit_concept_name %in% c("Other ambulatory visit",
                                                                  "Administrative Visit",
                                                                  "Interactive Telemedicine Service",
                                                                  "Non-hospital institution Visit",
                                                                  "Other") ~ "admin_other_telemed"
    )) %>% 
    group_by(person_id, visit_type_rough, visit_start_date) %>% 
    slice_min(visit_start_date, with_ties = FALSE) %>% 
    ungroup() %>% 
    group_by(person_id, visit_type_rough) %>% 
    summarise(n_visits_type = n()) %>% 
    ungroup() %>% 
    compute_new()
  
  visit_type_quantiles <-
    visit_lookback_counts %>% 
    group_by(visit_type_rough) %>% 
    summarise(quant20 = quantile(n_visits_type, probs=0.2),
              quant40 = quantile(n_visits_type, probs = .4), 
              quant60 = quantile(n_visits_type, probs = .6), 
              quant80 = quantile(n_visits_type, probs = 1)) %>% 
    ungroup() %>% 
    compute_new()
  
  
  cohort_with_util_flags <-
    visit_lookback_counts %>% 
    left_join(visit_type_quantiles, by="visit_type_rough") %>% 
    mutate(personal_visit_category = case_when(
      n_visits_type >= quant80 ~ paste("Num visits >=", quant80), 
      n_visits_type >= quant60 ~ paste("Num visits >=", quant60), 
      n_visits_type >= quant40 ~ paste("Num visits >=", quant40), 
      n_visits_type >= quant20 ~ paste("Num visits >=", quant20),
      n_visits_type < quant20 ~ paste("Num visits <", quant20)
    )) %>% 
    pivot_wider(id_cols = person_id, 
                names_from = visit_type_rough, 
                names_prefix = "visit_type_", 
                values_from = personal_visit_category,
                values_fill ="No visits")
  
  cohort %>% 
    left_join(cohort_with_util_flags, by="person_id") %>% 
    return()
  
}



#' Title
#' Quick summary function to report out the counts as fractions of columns the user specifies
#'
#' @param cohort_tbl 
#' @param group_column_1 
#' @param group_column_2 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' sum_fun <- function(.data, group, sum_var) {

summarise_as_fraction <- function(.data, group, distinct_key) {
  .data %>% 
    group_by(across({{ group }})) %>% 
    summarize("sum_{{distinct_key}}" := n_distinct({{ distinct_key }})) %>% 
    ungroup() %>% 
    mutate(total_sum = sum(sum_person_id)) %>% 
    mutate(prop = sum_person_id/total_sum) %>% 
    return()
}
  












