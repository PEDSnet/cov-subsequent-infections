
#' **Primary use function**
#'# Return a flag if one or more visits in the n years prior
#' Needs cohort to have column index_date (cohort entry date)
#' @param cohort 
#' @param visit_tbl Normally visit_ocurrence table
#' @param n_years Number of years to filer for a visit between index_date and n years ago
#'
#' @return Cohort table with new column `visit_within_n_years` either 1 or 0
#' @export
#'
#' @examples
flag_visits_prior <- function(cohort, visit_tbl, n_years) {
  
  prior_visits <- 
    visit_tbl %>% 
    select(person_id, visit_start_date) %>% 
    inner_join(cohort %>% 
                 select(person_id, index_date), by="person_id") %>% 
    group_by(person_id) %>% 
    filter(visit_start_date < index_date) %>% 
    slice_max(visit_start_date, with_ties = FALSE) %>% 
    mutate(visit_within_n_years = index_date - visit_start_date <= (365*n_years)) %>% 
    compute_new(indexes=c("person_id"))
  
  cohort %>% 
    left_join(prior_visits, by=c("person_id", "index_date")) %>% 
    select(-visit_start_date) %>% 
    mutate(visit_within_n_years = ifelse(is.na(visit_within_n_years), 0, 1)) %>% 
    return()
}

#' **Primary use function**
#' Return a flag if at least one visit MORE than 30 to 180 days after first index event (first covid positive)
#'
#' @param cohort 
#' @param visit_tbl Usually should be the visit_occurrence table
#'
#' @return Cohort table with new boolean columns stating whether a 
#' patient has had a visit at least n days following their index date 
#' @export
#'
#' @examples
flag_visit_later <- function(cohort, visit_tbl) {
  
  post_visits <- visit_tbl %>% 
    select(person_id, visit_start_date) %>% 
    inner_join(cohort %>% 
                 select(person_id, index_date), by="person_id") %>% 
    filter(visit_start_date > index_date) %>% 
    group_by(person_id) %>% 
    slice_max(visit_start_date, with_ties = FALSE) %>% # Person's most recent visit
    mutate(visit_atleast_30_days_post = (visit_start_date - index_date) >= 30) %>% ##Needs to be more than 30 days after the index date...
    mutate(visit_atleast_60_days_post = (visit_start_date - index_date) >= 60) %>% ##Needs to be more than 60 days after the index date...
    mutate(visit_atleast_90_days_post = (visit_start_date - index_date) >= 90) %>% ##Needs to be more than 90 days after the index date...
    mutate(visit_atleast_120_days_post = (visit_start_date - index_date) >= 120) %>% ##Needs to be more than 90 days after the index date...
    mutate(visit_atleast_150_days_post = (visit_start_date - index_date) >= 150) %>% ##Needs to be more than 90 days after the index date...
    mutate(visit_atleast_180_days_post = (visit_start_date - index_date) >= 180) %>% ##Needs to be more than 90 days after the index date...
    compute_new(indexes=c("person_id"))
  
  cohort %>% 
    left_join(post_visits, by=c("person_id", "index_date")) %>% 
    select(-visit_start_date) %>% 
    return()
}

flag_visit_follow_up <- function(cohort, visit_tbl) {
  
  folup_visits <- visit_tbl %>% 
    select(person_id, visit_start_date) %>% 
    inner_join(cohort %>% 
                 select(person_id, index_date), by="person_id") %>% 
    filter(visit_start_date >= as.Date(index_date + days(30)), visit_start_date < as.Date(index_date + days(180))) %>% 
    group_by(person_id) %>% 
    slice_max(visit_start_date, with_ties = FALSE) %>% # Person's most recent visit
    mutate(has_visit_in_followup_period = 1) %>% 
    # mutate(visit_atleast_30_days_post = (visit_start_date - index_date) >= 30) %>% ##Needs to be more than 30 days after the index date...
    # mutate(visit_atleast_60_days_post = (visit_start_date - index_date) >= 60) %>% ##Needs to be more than 60 days after the index date...
    # mutate(visit_atleast_90_days_post = (visit_start_date - index_date) >= 90) %>% ##Needs to be more than 90 days after the index date...
    # mutate(visit_atleast_120_days_post = (visit_start_date - index_date) >= 120) %>% ##Needs to be more than 90 days after the index date...
    # mutate(visit_atleast_150_days_post = (visit_start_date - index_date) >= 150) %>% ##Needs to be more than 90 days after the index date...
    # mutate(visit_atleast_180_days_post = (visit_start_date - index_date) >= 180) %>% ##Needs to be more than 90 days after the index date...
    compute_new(indexes=c("person_id"))
  
  cohort %>% 
    left_join(folup_visits, by=c("person_id", "index_date")) %>% 
    select(-visit_start_date) %>% 
    mutate(has_visit_in_followup_period = ifelse(is.na(has_visit_in_followup_period), 0, has_visit_in_followup_period)) %>% 
    return()
}
