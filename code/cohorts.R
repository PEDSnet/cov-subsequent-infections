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
#'
apply_inclusion_flags_rsv_study <- function(cohort_tbl, ce_start_date, ce_end_date, cohort_label) {
  eligible_pats <-
    cohort_tbl %>% 
    mutate(index_date = ce_date) %>% 
    flag_visits_prior(cdm_tbl("visit_occurrence"), n_years = 1.5) %>% 
    flag_visit_follow_up(cdm_tbl("visit_occurrence")) %>% 
    select(-index_date) %>% 
    compute_new(indexes=c("person_id")) %>% 
    flag_study_eligiblity(ce_start_date, ce_end_date) %>% 
    compute_new(indexes=c("person_id"))
  
  eligible_pats %>% 
    return()
}

pats_with_visit_under_age <- function(person_tbl, visit_tbl, ce_start_date, ce_end_date, age_limit_years_rough) {
  # Gets all patients under 6 during the whole study entry period
  person_tbl %>% 
    select(person_id, birth_date) %>% 
    filter(as.numeric(ce_end_date - birth_date)/365.25 <= 6) %>% 
    inner_join(visit_tbl %>% 
                 filter(visit_start_date >= ce_start_date, visit_start_date < ce_end_date) %>% 
                 select(person_id, visit_start_date), by="person_id") %>% 
    group_by(person_id) %>% 
    slice_min(visit_start_date, with_ties = FALSE) %>% 
    ungroup() %>% 
    select(person_id, birth_date) %>% 
    return()
}

get_influenza_evidence <- function(cohort_tbl, ce_start_date, ce_end_date) {
  ## Return the cohort with a flu evidence within the study period
  ## columns: person_id, flu_date
  earliest_date <- as.Date(ce_start_date - months(12))
  latest_date <- as.Date(ce_end_date + months(6))
  influenza_positives <-
    cohort_tbl %>%
    get_flu_evidence(min_date = earliest_date, max_date = latest_date)
  
  cohort_tbl %>% 
    left_join(influenza_positives %>% 
                select(person_id, 
                       flu_date = earliest_flu_evidence,
                       flu_dx_date,
                       earliest_flu_test
                       ),
              by="person_id") %>% 
    return()
}

get_influenza_evidence_test_only <- function(cohort_tbl, ce_start_date, ce_end_date) {
  ## Return the cohort with a flu evidence within the study period
  ## columns: person_id, flu_date
  earliest_date <- as.Date(ce_start_date - months(12))
  latest_date <- as.Date(ce_end_date + months(6))
  influenza_positives <-
    cohort_tbl %>%
    get_flu_evidence(min_date = earliest_date, max_date = latest_date, test_only = TRUE)
  
  cohort_tbl %>% 
    left_join(influenza_positives %>% 
                select(person_id, 
                       flu_date = earliest_flu_evidence,
                       earliest_flu_test,
                       flu_lab_concept_id,
                       flu_value_concept_id
                ),
              by="person_id") %>% 
    return()
}

get_other_resp_evidence <- function(cohort_tbl, ce_start_date, ce_end_date) {
  
  earliest_date <- as.Date(ce_start_date - months(12))
  latest_date <- as.Date(ce_end_date + months(6))
  
  respiratory_positives <-
    cohort_tbl %>%
    get_respiratory_evidence(min_date = earliest_date, max_date = latest_date)
  
  cohort_tbl %>% 
    left_join(respiratory_positives %>% 
                select(person_id, 
                       resp_date = condition_start_date,
                       resp_concept = concept_name,
                       resp_concept_id = concept_id
                ),
              by="person_id") %>% 
    return()
}

get_covid_evidence <- function(cohort_tbl, odr_tbl, ce_start_date, ce_end_date) {
  ## First, filter ODR table to the max post-acute and pre-acute period as necessary for the exclusion criteria 
  earliest_date <- as.Date(ce_start_date - months(12))
  latest_date <- as.Date(ce_end_date + months(6))
  odr_padded <- odr_tbl %>% 
    filter(as.Date(observation_date) >= as.Date("2020-03-01"), # TODO actually maybe don't padd this, because patients need to be washed out of having any prior covid
           as.Date(observation_date) < latest_date) %>% 
    compute_new()
  
  covid_positives <- 
    cohort_tbl %>% 
    flag_covid_positives(odr_tbl = odr_padded)
  
  cohort_tbl %>% 
    left_join(covid_positives %>% 
                select(person_id, covid_date = index_date,
                       covid_index_date_imputed = index_date_imputed,
                       cov_test_date = test_date, covid_dx_date, pasc_dx_date),
              by="person_id") %>% 
    return()
  
}

get_covid_evidence_test_only <- function(cohort_tbl, odr_tbl, ce_start_date, ce_end_date) {
  ## First, filter ODR table to the max post-acute and pre-acute period as necessary for the exclusion criteria 
  earliest_date <- as.Date(ce_start_date - months(12))
  latest_date <- as.Date(ce_end_date + months(6))
  odr_padded <- odr_tbl %>% 
    filter(as.Date(observation_date) >= earliest_date, 
           as.Date(observation_date) < latest_date) %>% 
    compute_new()
  
  covid_positives <- 
    cohort_tbl %>% 
    flag_covid_positives_test_only(odr_tbl = odr_padded) %>% 
    mutate(covid_obs_concept_id = observation_concept_id, covid_value_concept_id = value_as_concept_id)
  
  cohort_tbl %>% 
    left_join(covid_positives %>% 
                select(person_id,
                       covid_date = test_date,
                       covid_obs_concept_id,
                       covid_value_concept_id),
              by="person_id") %>% 
    return()
  
}

### This will be the primary use function
build_comparison_cohorts_rsv_study <- function(cohort_tbl, odr_tbl, ce_start_date, ce_end_date) {
  # Logic in here to decide what to do with kids who had both within the study period, or who had one or the other or none
  relevant_infections <-
    sample_pats_under_6 %>% #cohort_tbl %>% 
    get_covid_evidence(odr_tbl, ce_start_date, ce_end_date) %>% 
    get_influenza_evidence(ce_start_date, ce_end_date) %>% 
    get_other_resp_evidence(ce_start_date, ce_end_date) %>% ## shouldn't there be data here on pre- and post- cohort entry period infections?
    compute_new(indexes=c("person_id"))
  
  ## Then relevant infections should have the earliest infection for each person, in each category of comparison cohort
  ## Next step, apply exclusion rules
  ## Apply rules:
  ## If one of the other categories has happened 12-6 months before index period, count as variable
  ## IF one of the other categories has happened within 6 months of index, exclude
  ## If one of the other categories has happened within 6 months following index, exclude
  ## If other category is after post-acute period, ok
  ## new columns: comparison_cohort, prior_infection_type, exclude (1/0), exclude_reason
  ## classify phases instead (will make logic easier later)
  # relevant_infections <- data.frame(person_id = c("a", "b", "c", "d", "e", "f"), 
  #                                   covid_date = c(as.Date("2022-04-01"), as.Date("2022-04-05"), NA, as.Date("2021-09-01"), as.Date("2022-04-01"), NA), 
  #                                   resp_date = c(as.Date("2022-06-01"), NA, as.Date("2022-07-01"), as.Date("2022-03-02"), as.Date("2022-04-01"), NA), 
  #                                   flu_date = c(as.Date("2024-01-01"), as.Date("2022-08-01"), NA, NA, NA, NA),
  #                                   birth_date = c(as.Date("2018-03-01"), as.Date("2009-01-01"), as.Date("2015-01-01"),as.Date("2016-01-01"),as.Date("2008-01-01"),as.Date("2021-01-01")))
  # 
  infections_phases_categorized <-
    relevant_infections %>% 
    pivot_longer(cols = c("covid_date", "flu_date", "resp_date"), names_to = "event", values_to = "event_date") %>% 
    mutate(event_in_ce_window = ifelse(event_date >= ce_start_date & event_date < ce_end_date, 1, 0)) %>% 
    select(person_id, cov_test_date, covid_dx_date, flu_dx_date, earliest_flu_test, event, event_date,
           event_in_ce_window) %>% 
    filter(!is.na(event_date)) %>% 
    mutate(lab_confirmed = case_when(
      event=="covid_date" & !is.na(cov_test_date) ~ 1,
      event=="flu_date" & !is.na(earliest_flu_test) ~ 1,
      TRUE ~ 0
    )) %>% 
    compute_new(indexes=c("person_id"))
  
  ## Now need to identify the earliest index event in the window
  ## Anyone with NO evidence of all 3 will be removed during this step because they have no minimum
  ## New step: earliest index date, test confirmed, and earliest index date, any (I guess it's just, what is the earliest thing, and if its test confirmed or not)
  earliest_index_event <-
    infections_phases_categorized %>% 
    filter(event_in_ce_window == 1) %>% 
    group_by(person_id) %>% 
    slice_min(event_date, with_ties = FALSE) %>% 
    ungroup() %>% 
    select(person_id, ce_date = event_date, ce_event = event, lab_confirmed) %>%  
    compute_new(indexes=c("person_id"))
    
  infections_phases_relative <-
    earliest_index_event %>% 
    left_join(relevant_infections, by="person_id") %>% 
    mutate(covid_phase = case_when(
      is.na(covid_date) ~ "none",
      ce_event == "covid_date" ~ "index",
      covid_date < as.Date(ce_date - months(6)) ~ "prior",
      covid_date < ce_date ~ "pre-acute",
      covid_date < as.Date(ce_date + months(1)) ~ "acute",
      covid_date < as.Date(ce_date + months(6)) ~ "post-acute",
      covid_date >= as.Date(ce_date + months(6)) ~ "future"
    )) %>% 
    mutate(flu_phase = case_when(
      is.na(flu_date) ~ "none",
      ce_event == "flu_date" ~ "index",
      flu_date < as.Date(ce_date - months(6)) ~ "prior",
      flu_date < ce_date ~ "pre-acute",
      flu_date < as.Date(ce_date + months(1)) ~ "acute",
      flu_date < as.Date(ce_date + months(6)) ~ "post-acute",
      flu_date >= as.Date(ce_date + months(6)) ~ "future"
    )) %>% 
    mutate(resp_phase = case_when(
      is.na(resp_date) ~ "none",
      ce_event == "resp_date" ~ "index",
      resp_date < as.Date(ce_date - months(6)) ~ "prior",
      resp_date < ce_date ~ "pre-acute",
      resp_date < as.Date(ce_date + months(1)) ~ "acute",
      resp_date < as.Date(ce_date + months(6)) ~ "post-acute",
      resp_date >= as.Date(ce_date + months(6)) ~ "future"
    )) 
  
  exclude_reasons <-
    infections_phases_relative %>% 
    mutate(exclude = case_when(
      covid_phase %in% c("prior","pre-acute", "acute", "post-acute") |
        flu_phase %in% c("pre-acute", "acute", "post-acute") ~ 1,
        resp_phase %in% c("pre-acute", "acute", "post-acute") ~ 2, # Changing this to it's ok if someone has a diagnosis of a respiratory condition
      TRUE ~ 0
    )) %>% 
    mutate(exclude_reason = case_when(
      exclude == 2 ~ "Has respiratory dx within pre-acute, acute, or post-acute period",
      exclude == 0 ~ "none",
      TRUE ~ paste0("Had COVID: ", covid_phase, ", FLU: ", flu_phase)
    ))
  
  exclude_reasons %>% 
    mutate(sub_cohort = case_when(
      ce_event == "covid_date" ~ "COVID",
      ce_event == "flu_date" ~ "Influenza",
      ce_event == "resp_date" ~ "Respiratory")) %>% 
    mutate(age_years_on_ce_date = as.numeric(ce_date-birth_date)/365.25) %>% 
    compute_new(indexes=c("person_id")) %>% 
    return()
  
}

build_comparison_cohorts_rsv_study_sensitivity <- function(cohort_tbl, odr_tbl, ce_start_date, ce_end_date) {
  # Logic in here to decide what to do with kids who had both within the study period, or who had one or the other or none
  # TODO change this to be the same function, regular code, with just a flag at the end for test-only or not, so the cohort can actually just be filtered for evidence
  relevant_infections <-
    cohort_tbl %>% 
    get_covid_evidence_test_only(odr_tbl, ce_start_date, ce_end_date) %>% 
    get_influenza_evidence_test_only(ce_start_date, ce_end_date) %>% 
    get_other_resp_evidence(ce_start_date, ce_end_date) %>% ## shouldn't there be data here on pre- and post- cohort entry period infections?
    compute_new(indexes=c("person_id"))
  
  ## Then relevant infections should have the earliest infection for each person, in each category of comparison cohort
  ## Next step, apply exclusion rules
  ## Apply rules:
  ## If one of the other categories has happened 12-6 months before index period, count as variable
  ## IF one of the other categories has happened within 6 months of index, exclude
  ## If one of the other categories has happened within 6 months following index, exclude
  ## If other category is after post-acute period, ok
  ## new columns: comparison_cohort, prior_infection_type, exclude (1/0), exclude_reason
  ## classify phases instead (will make logic easier later)
  # relevant_infections <- data.frame(person_id = c("a", "b", "c", "d", "e", "f"), 
  #                                   covid_date = c(as.Date("2022-04-01"), as.Date("2022-04-05"), NA, as.Date("2021-09-01"), as.Date("2022-04-01"), NA), 
  #                                   resp_date = c(as.Date("2022-06-01"), NA, as.Date("2022-07-01"), as.Date("2022-03-02"), as.Date("2022-04-01"), NA), 
  #                                   flu_date = c(as.Date("2024-01-01"), as.Date("2022-08-01"), NA, NA, NA, NA),
  #                                   birth_date = c(as.Date("2018-03-01"), as.Date("2009-01-01"), as.Date("2015-01-01"),as.Date("2016-01-01"),as.Date("2008-01-01"),as.Date("2021-01-01")))
  # 
  infections_phases_categorized <-
    relevant_infections %>% 
    pivot_longer(cols = c("covid_date", "flu_date", "resp_date"), names_to = "event", values_to = "event_date") %>% 
    mutate(event_in_ce_window = ifelse(event_date >= ce_start_date & event_date < ce_end_date, 1, 0)) %>% 
    compute_new(indexes=c("person_id"))
  
  ## Now need to identify the earliest index event in the window
  ## Anyone with NO evidence of all 3 will be removed during this step because they have no minimum
  earliest_index_event <-
    infections_phases_categorized %>% 
    filter(event_in_ce_window == 1) %>% 
    group_by(person_id) %>% 
    slice_min(event_date, with_ties = FALSE) %>% 
    ungroup() %>% 
    select(person_id, ce_date = event_date, ce_event = event) %>%  
    compute_new(indexes=c("person_id"))
  
  infections_phases_relative <-
    earliest_index_event %>% 
    left_join(relevant_infections, by="person_id") %>% 
    mutate(covid_phase = case_when(
      is.na(covid_date) ~ "none",
      ce_event == "covid_date" ~ "index",
      covid_date < as.Date(ce_date - months(6)) ~ "prior",
      covid_date < ce_date ~ "pre-acute",
      covid_date < as.Date(ce_date + months(1)) ~ "acute",
      covid_date < as.Date(ce_date + months(6)) ~ "post-acute",
      covid_date >= as.Date(ce_date + months(6)) ~ "future"
    )) %>% 
    mutate(flu_phase = case_when(
      is.na(flu_date) ~ "none",
      ce_event == "flu_date" ~ "index",
      flu_date < as.Date(ce_date - months(6)) ~ "prior",
      flu_date < ce_date ~ "pre-acute",
      flu_date < as.Date(ce_date + months(1)) ~ "acute",
      flu_date < as.Date(ce_date + months(6)) ~ "post-acute",
      flu_date >= as.Date(ce_date + months(6)) ~ "future"
    )) %>% 
    mutate(resp_phase = case_when(
      is.na(resp_date) ~ "none",
      ce_event == "resp_date" ~ "index",
      resp_date < as.Date(ce_date - months(6)) ~ "prior",
      resp_date < ce_date ~ "pre-acute",
      resp_date < as.Date(ce_date + months(1)) ~ "acute",
      resp_date < as.Date(ce_date + months(6)) ~ "post-acute",
      resp_date >= as.Date(ce_date + months(6)) ~ "future"
    )) 
  
  exclude_reasons <-
    infections_phases_relative %>% 
    mutate(exclude = case_when(
      covid_phase %in% c("pre-acute", "acute", "post-acute") |
        flu_phase %in% c("pre-acute", "acute", "post-acute") |
        resp_phase %in% c("pre-acute", "acute", "post-acute") ~ 1,
      TRUE ~ 0
    )) %>% 
    mutate(exclude_reason = case_when(
      exclude == 0 ~ "none",
      TRUE ~ paste0("Had COVID: ", covid_phase, ", FLU: ", flu_phase, ", RESP: ", resp_phase)
    ))
  
  exclude_reasons %>% 
    mutate(sub_cohort = case_when(
      ce_event == "covid_date" ~ "COVID",
      ce_event == "flu_date" ~ "Influenza",
      ce_event == "resp_date" ~ "Respiratory")) %>% 
    mutate(age_years_on_ce_date = as.numeric(ce_date-birth_date)/365.25) %>% 
    compute_new(indexes=c("person_id")) %>% 
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


### Performing the washouts
## Getting data pertaining to the washouts
#MISC or PASC: if observation_concept_id==2000001527, value_as_concept_id==703578, or 2000001520  and patient does NOT have covid, then washout
# Serology test, observation_concept_id==2000001528, value_as_concept_id %in% c(2000001526, 9191)
#' Title
#'
#' @param cohort 
#' @param odr 
#'
#' @return
#' @export
#'
#' @examples
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

cohort_demo <- function(cohort_tbl) {
  eligible_pats_demo <-
    cohort_tbl %>% 
    left_join(cdm_tbl("person") %>% 
                select(person_id, gender_concept_id,
                       race_concept_id, ethnicity_concept_id), by="person_id") %>% 
    mutate(mock_medical_complexity = case_when(person_id %% 17 == 1 ~ "complex chronic",
                                               person_id %% 8==1 ~ "chronic, non-complex",
                                               TRUE ~ "non complex or chronic")) %>% 
    mutate(sex_cat = case_when(gender_concept_id == 8507L ~ 'Male',
                               gender_concept_id == 8532L ~ 'Female',
                               TRUE ~ 'Other/unknown'),
           race_eth_cat = case_when(ethnicity_concept_id == 38003563L ~ 'Hispanic',
                                    race_concept_id == 8516L ~ 'NH_Black/AA',
                                    race_concept_id %in% c(8515L, 8557L) ~
                                      'NH_Asian/PI',
                                    #                               race_concept_id == 8657L ~ 'Native American',
                                    race_concept_id == 8527L ~ 'NH_White',
                                    race_concept_id == 44814659L ~ 'NH_Multiple',
                                    TRUE ~ 'Other/Unknown'),
           age_cat_at_ce = case_when(age_years_on_ce_date < 0.5 ~ "<6 months",
                                     age_years_on_ce_date < 1 ~ "6m-1y old",
                                     age_years_on_ce_date < 2 ~ "1-2y old",
                                     age_years_on_ce_date < 3 ~ "2-3y old",
                                     age_years_on_ce_date < 4 ~ "3-4y old",
                                     TRUE ~ "4-5y old")) %>% 
    mutate(cohort_entry_week = floor_date(ce_date, unit="week")) %>% 
    select(-gender_concept_id, -race_concept_id, -ethnicity_concept_id) %>% 
    flag_utilization_level(cdm_tbl("visit_occurrence"), lookback_days = 365) %>% 
    # left_join(results_tbl(paste0("cohort", cohort_1_label, "pmca")), by="person_id") %>% 
    compute_new(indexes=c("person_id"))
  
  # eligible_pats_demo %>% 
  #   output_tbl(cohort_output_tbl_name, indexes=c("person_id"))
  # 
  eligible_pats_demo %>% 
    return()
}


generate_output_cohort_demo <- function(ce_start_date, ce_end_date, cohort_tbl_name, cohort_output_tbl_name) {
  eligible_pats_demo <-
    results_tbl(cohort_tbl_name) %>% 
    left_join(cdm_tbl("person") %>% 
                select(person_id, gender_concept_id,
                       race_concept_id, ethnicity_concept_id), by="person_id") %>% 
    mutate(mock_medical_complexity = case_when(person_id %% 17 == 1 ~ "complex chronic",
                                               person_id %% 8==1 ~ "chronic, non-complex",
                                               TRUE ~ "non complex or chronic")) %>% 
    mutate(sex_cat = case_when(gender_concept_id == 8507L ~ 'Male',
                               gender_concept_id == 8532L ~ 'Female',
                               TRUE ~ 'Other/unknown'),
           race_eth_cat = case_when(ethnicity_concept_id == 38003563L ~ 'Hispanic',
                                    race_concept_id == 8516L ~ 'NH_Black/AA',
                                    race_concept_id %in% c(8515L, 8557L) ~
                                      'NH_Asian/PI',
                                    #                               race_concept_id == 8657L ~ 'Native American',
                                    race_concept_id == 8527L ~ 'NH_White',
                                    race_concept_id == 44814659L ~ 'NH_Multiple',
                                    TRUE ~ 'Other/Unknown'),
           age_cat_at_ce = case_when(age_years_on_ce_date < 0.5 ~ "<6 months",
                                     age_years_on_ce_date < 1 ~ "6m-1y old",
                                     age_years_on_ce_date < 2 ~ "1-2y old",
                                     age_years_on_ce_date < 3 ~ "2-3y old",
                                     age_years_on_ce_date < 4 ~ "3-4y old",
                                     TRUE ~ "4-5y old")) %>% 
    mutate(cohort_entry_week = floor_date(ce_date, unit="week")) %>% 
    select(-gender_concept_id, -race_concept_id, -ethnicity_concept_id) %>% 
    flag_utilization_level(cdm_tbl("visit_occurrence"), lookback_days = 365) %>% 
    # left_join(results_tbl(paste0("cohort", cohort_1_label, "pmca")), by="person_id") %>% 
    compute_new(indexes=c("person_id"))
  
  eligible_pats_demo %>% 
    output_tbl(cohort_output_tbl_name, indexes=c("person_id"))
  
  eligible_pats_demo %>% 
    return()
}

generate_output_outcome_lists <- function(cohort, outcome, outcome_start_date, outcome_end_date, output_tbl_name, test_only = FALSE) {
  if (outcome == "rsv") {
    outcomes <- cohort %>% 
      flag_rsv_outcome(outcome_start_date, outcome_end_date, test_only) %>% 
      compute_new(indexes=c("person_id"))
  } else if (outcome == "general") {
    outcomes <- cohort %>% 
      flag_outcome_infections(outcome_start_date,
                              outcome_end_date) %>% 
      compute_new(indexes=c("person_id"))
  } else if (outcome == "respiratory") {
    outcomes <- cohort %>% 
      flag_outcome_resp_infections(outcome_start_date,
                                   outcome_end_date) %>% 
      compute_new(indexes=c("person_id"))
  } else {
    print(paste0("Hm, did not have a method available for ", outcome," , available outcomes are rsv, general, and respiratory."))
  }
  
  outcomes %>% 
    output_tbl(output_tbl_name)
}

### Function to compute trajectories for patient outcomes
# compute_patient_trajectories(cohort) {
#   
# }


plot_rollup_summary <- function(rolled_up_dx_tbl, 
                                n_sample_size = 50, 
                                person_labels = FALSE, show_intermed_events = TRUE) {
  if (show_intermed_events) {
    data <- rolled_up_dx_tbl
  } else {
    data <- rolled_up_dx_tbl %>% 
      filter(record_flag != "intermediate_follow_up_dx")
  }
  
  rollup_plot <-
    data %>% 
    group_by(person_id) %>% 
    mutate(max_day = max(days_from_ce)) %>% 
    mutate(viz_id=cur_group_id()) %>% 
    ungroup() %>%
    mutate(viz_id = fct_reorder(factor(viz_id), max_day)) %>% 
    group_by(person_id, event) %>% 
    mutate(person_event = paste0(person_id, "_", event)) %>% 
    # ungroup() %>% 
    ggplot() +
    geom_vline(xintercept = -365, linetype="dashed", color="grey") +
    geom_vline(xintercept = -182.5, linetype="dashed", color="grey") +
    geom_vline(xintercept = 30, linetype="dashed", color="orange") +
    geom_vline(xintercept = 180, linetype="dashed", color="orange") +
    geom_vline(xintercept = 0, linetype="solid", size=0.5) +
    geom_line(aes(x=days_from_ce, y = viz_id, group=person_event),
              linewidth = 1, alpha = 0.9) +
    geom_line(aes(x=days_from_ce, y = viz_id, group=viz_id),
              linewidth = 0.25, alpha = 0.2) +
    geom_point(aes(x=days_from_ce, y = viz_id, shape = record_flag, col=record_flag),
               stroke = 1.5,  alpha = .5) +
    theme_bw() +
    scale_color_manual(values = wes_palette("FantasticFox1", n=8, type = "continuous")) +
    labs(x="Days from cohort entry date", y= "Patient ID")
  
  if (person_labels) {
    rollup_plot <- rollup_plot + 
      geom_label(aes(x=-10, y = viz_id, group=viz_id, label=person_id))
  }
  
  rollup_plot %>% 
    return()
}


generate_outcome_dataset <- function(cohort,
                                     outcome_type,
                                     respiratory_tbl,
                                     general_tbl,
                                     rsv_tbl,
                                     output_tbl_name) {
  if (outcome_type=="presence") {
    ## Then roll up the data into earliest event of that type
    earliest_resp_infection <-
      respiratory_tbl %>% 
      left_join(cohort, by="person_id") %>% 
      mutate(days_til_resp_outcome = as.numeric(condition_start_date - ce_date)) %>% 
      filter(days_til_resp_outcome > 7) %>% 
      group_by(person_id) %>% 
      slice_min(condition_start_date, with_ties = FALSE) %>% 
      select(person_id, earliest_resp_outcome = condition_start_date, resp_concept = concept_name) %>% 
      ungroup() %>% 
      compute_new()
    
    cohort_with_resp_outcomes <-
      cohort %>% 
      left_join(earliest_resp_infection, by="person_id") %>% 
      mutate(days_til_earliest_resp_outcome = as.numeric(earliest_resp_outcome - ce_date)) %>% 
      compute_new()

    earliest_general_infection <-
      general_tbl %>% 
      left_join(cohort, by="person_id") %>% 
      mutate(days_til_any_outcome = as.numeric(condition_start_date - ce_date)) %>% 
      filter(days_til_any_outcome > 7) %>% 
      group_by(person_id) %>% 
      slice_min(condition_start_date, with_ties = FALSE) %>% 
      select(person_id, earliest_general_outcome = condition_start_date, general_concept = concept_name) %>% 
      ungroup() %>% 
      compute_new()
    
    cohort_with_general_outcomes <-
      cohort_with_resp_outcomes %>% 
      left_join(earliest_general_infection, by="person_id") %>% 
      mutate(days_til_earliest_general_outcome = as.numeric(earliest_general_outcome - ce_date)) %>% 
      compute_new()
    
    cohort_with_all_outcomes <-
      cohort_with_general_outcomes %>% 
      left_join(rsv_tbl %>% 
                  select(person_id, rsv_evidence_date), by = "person_id") %>% 
      mutate(days_til_earliest_rsv_outcome = as.numeric(rsv_evidence_date - ce_date)) %>% 
      mutate(has_postacute_resp_outcome = case_when(
        !is.na(resp_concept) & days_til_earliest_resp_outcome >= 30 & days_til_earliest_resp_outcome < 180 ~ 1,
        TRUE ~ 0
      )) %>% 
      mutate(has_postacute_general_outcome = case_when(
        !is.na(general_concept) & days_til_earliest_general_outcome >= 30 & days_til_earliest_general_outcome < 180 ~ 1,
        TRUE ~ 0
      )) %>% 
      mutate(has_postacute_rsv_outcome = case_when(
        !is.na(rsv_evidence_date) & days_til_earliest_rsv_outcome >= 30 & days_til_earliest_rsv_outcome < 180 ~ 1,
        TRUE ~ 0
      )) %>% 
      compute_new()
    
    cohort_with_all_outcomes %>% 
      output_tbl(output_tbl_name)
    
    cohort_with_all_outcomes %>% 
      return()
    
  } else if (outcome_type == "count") {
    ## Generate total number of specific events per person
    ### this should include deduplication logic, can go into sub-functions
  }
  
}
  

get_insurance_class <- function(cohort) {
  insurance_plans <- cohort %>% 
    left_join(cdm_tbl("visit_occurrence") %>% 
                select(visit_occurrence_id, person_id, visit_start_date), by=c("person_id", "ce_date"="visit_start_date")) %>%
    compute_new() %>% 
    left_join(cdm_tbl("visit_payer") %>% 
                select(visit_occurrence_id, plan_class), by="visit_occurrence_id") %>% 
    mutate(earliest_date = ce_date) %>% 
    distinct(person_id, plan_class, earliest_date) %>% 
    ungroup() %>% 
    compute_new()
  
  insurance_plans_labeled <-
    insurance_plans %>% 
    mutate(plan_class_index=case_when(
      plan_class %in% c("Other public", "Medicaid/sCHIP", "Medicare")~2,
      plan_class %in% c("Private/Commercial")~1,
      plan_class %in% c("Self-pay", "Other/Unknown")~0
    )) %>%
    group_by(person_id) %>%
    slice_max(plan_class_index, with_ties=FALSE) %>%
    mutate(insurance_class=case_when(
      plan_class_index==2~"public",
      plan_class_index==1~"private",
      plan_class_index==0~"other/unknown"
    )) %>% 
    compute_new(indexes=c("person_id"))
  
  cohort %>% 
    left_join(insurance_plans_labeled %>% 
                select(person_id, insurance_class), by="person_id") %>% 
    return()
  
}












