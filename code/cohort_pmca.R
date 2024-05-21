#' Function from pasc25_hazard_ratio / code / cohort_pmca
#' produce PMCA table with 3 year lookback
#'
#' @param cohort cohort of patients with patients and observation_date
#' @param pmca_xwalk codeset that has flags for body systems and whether or not progressive
#' @param condition_tbl formatting of the condition occurrence table
#' @param visit_tbl formatting of the visit occurrence table
#'
#' @return table that has conditions and flags for body systems
#' columns:
#' person_id | observation_date | result_derivation | condition_concept_id | condition_concept_name |
#' condition_source_concept_id | condition_source_value | visit_occurrence_id | condition_start_state |
#' description | body_system | progressive | visit_concept_id | site
#'

produce_pmca_lookup <- function (cohort,
                                 pmca_xwalk=load_codeset('pmca_icd10'),
                                 condition_tbl=cdm_tbl('condition_occurrence'),
                                 visit_tbl=cdm_tbl('visit_occurrence')) {
  only_pmca_conds <-
    cohort %>%
    inner_join(
      select(condition_tbl,
             person_id,
             condition_concept_id,
             condition_concept_name,
             condition_source_concept_id,
             condition_source_concept_name,
             visit_occurrence_id,
             condition_start_date),
      by=c('person_id')
    ) %>% filter(observation_date - condition_start_date <= 1096L &
                   observation_date - condition_start_date >= 0L) %>%
    
    inner_join(
      pmca_xwalk,
      by=c('condition_source_concept_id'='concept_id')
    ) %>%
    inner_join(
      select(
        visit_tbl,
        visit_occurrence_id,
        visit_concept_id
      )
    ) %>% filter(visit_concept_id %in% c(9202L,9201L,9203L,
                                         2000000048L,2000000088L,581399L)) %>%
    distinct() 
  
}



#' Function from pasc25_hazard_ratio / code / cohort_pmca
#' compute patient, body system with visit number, and flags for malignancy and progressive
#'
#' @param pmca_lookup_tbl pmca table output from `produce_pmca_lookup`;
#' must contain `observation_date` and `body_system` and `condition_start_date` and flagged conditions
#'
#' @return computes information to be able to apply algorithms; groups by body system and counts visits,
#' with flags for progressive or malignancy for patients
#'
#' person_id | observation_date | result_derivation | body_system |
#' yr_1 | yr_2 | yr_3 | total_visits | progressive | malignancy
#'

compute_pmca_summary <- function(pmca_lookup_tbl) {
  
  add_year <-
    pmca_lookup_tbl %>%
    filter(
      ! body_system == 'malignancy'
    ) %>% mutate(
      flag_yr = case_when(
        condition_start_date < observation_date &
          condition_start_date >= sql("(observation_date - interval '1 year')::date") ~ 'yr_1',
        condition_start_date < sql("(observation_date - interval '1 year')::date") &
          condition_start_date >= sql("(observation_date - interval '2 years')::date") ~ 'yr_2',
        condition_start_date < sql("(observation_date - interval '2 years')::date") &
          condition_start_date >= sql("(observation_date - interval '3 years')::date") ~ 'yr_3',
        TRUE ~ 'no_yr'
      )) %>%
    group_by(
      person_id,
      observation_date,
      #result_derivation,
      body_system,
      flag_yr
    ) %>% summarise(
      visit_yr_ct=as.integer(n_distinct(condition_start_date))
    ) %>% pivot_wider(names_from = flag_yr,
                      values_from = visit_yr_ct,
                      values_fill = 0L) %>%
    #collect() %>%
    #mutate(total_visits=sum(c_across(starts_with("yr_")), na.rm=T)) %>%
    mutate(total_visits = yr_1 + yr_2 + yr_3) %>%
    ungroup() %>% compute_new(temporary=TRUE,
                              indexes=list('person_id'))
  
  progressive_malignant_pts <-
    pmca_lookup_tbl %>%
    filter(
      progressive == 'yes' |
        body_system == 'malignancy'
    ) %>% mutate(malignancy =
                   case_when(
                     body_system == 'malignancy' ~ 'yes',
                     TRUE ~ 'no'
                   )) %>%
    filter(malignancy=='yes' | progressive == 'yes') %>%
    select(person_id,
           progressive,
           malignancy) %>% distinct %>%
    compute_new(temporary=TRUE,
                indexes=list('person_id'))
  
  all_pts <-
    dplyr::union(
      add_year %>% ungroup() %>% select(person_id),
      progressive_malignant_pts %>% select(person_id)
    )
  
  all_pts %>%
    left_join(add_year) %>%
    left_join(progressive_malignant_pts) %>%
    mutate(
      progressive=case_when(
        is.na(progressive) ~ 'no',
        TRUE ~ progressive
      ),
      malignancy=case_when(
        is.na(malignancy) ~ 'no',
        TRUE ~ malignancy
      )
    ) %>% select(
      person_id,observation_date,#result_derivation,
      body_system,yr_1,yr_2,yr_3,total_visits,progressive,malignancy
    )
  
}



#' compute classification algorithm for *most conservative*: 
#' *complex chronic* is defined as patients with at least 1 visit for two body systems for all three years OR progressive OR malignant
#' *chronic* is defined as having at least 1 visit all three years for just one body system
#' 
#' @param pmca_lookup_tbl output from `produce_pmca_lookup`
#' @param cohort_tbl cohort table with person_id, observation_date and site
#' 
#' @return table that has patients in the *most conservative* category with the following columns:
#' person_id | observation_date | body_system | yr_1 | yr_2 | yr_3 | total_visits |
#' progressive | malignancy | complex_chronic | chronic | non_complex_chronic
#'

compute_pmca_cats_cons <- function(pmca_summary_tbl,
                                   cohort_tbl) {
  
  gt_two_bs <- 
    pmca_summary_tbl %>%
    filter(
      yr_1 > 0 &
        yr_2 > 0 &
        yr_3 > 0
    ) %>%
    group_by(person_id,
             observation_date) %>%
    summarise(body_system_ct=n_distinct(body_system)) %>%
    filter(
      body_system_ct > 1
    ) %>% ungroup()
  
  prog_or_malig <-
    pmca_summary_tbl %>% 
    filter(progressive == 'yes' | malignancy == 'yes') 
  
  complex_pts <- 
    dplyr::union(
      gt_two_bs %>% select(person_id,observation_date),
      prog_or_malig %>% select(person_id,observation_date)
    ) %>% mutate(complex_chronic = 1L) %>% 
    compute_new(temporary=TRUE,
                indexes=list('person_id'))
  
  chronic_pts <- 
    pmca_summary_tbl %>%
    filter(
      yr_1 > 0 &
        yr_2 > 0 &
        yr_3 > 0
    ) %>% anti_join(complex_pts,
                    by=c('person_id','observation_date')) %>%
    distinct(person_id, observation_date) %>% mutate(chronic = 1L) %>%
    compute_new(temporary=TRUE,
                indexes=list('person_id'))
  
  pmca_summary_tbl %>%
    right_join(cohort_tbl,
               by=c('person_id','observation_date')) %>% 
    left_join(complex_pts) %>%
    left_join(chronic_pts) %>%
    mutate(complex_chronic=case_when(is.na(complex_chronic) ~ 0L,
                                     TRUE ~ complex_chronic),
           chronic=case_when(is.na(chronic) ~ 0L,
                             TRUE ~ chronic)) %>%
    mutate(non_complex_chronic = 
             case_when(complex_chronic == 1 | chronic == 1 ~ 0L,
                       TRUE ~ 1L)) %>% select(person_id,
                                              observation_date,
                                              site,
                                              complex_chronic,
                                              chronic,
                                              non_complex_chronic) %>% distinct()
}





#' @param pmca_summary pmca summary table output from compute_pmca_summary. Counts are unique by person_id, body_system, and malignancy/progressive.
#' @param cohort_tbl full cohort tbl
#' @param pmca_category pmca_category table output from compute_pmca_cat_cons function
#' @return table with person-level flags for progressive/malignant and number of body systems affected
#' 
pmca_all_flags <- function(pmca_summary_tbl,
                           cohort_tbl,
                           pmca_category) {
  
  progressive_malignancy_flag <- pmca_summary_tbl %>%
    mutate(progressive_num=case_when(progressive=='yes' ~ 1L,
                                     progressive=='no' ~ 0L),
           malignancy_num=case_when(malignancy=='yes' ~ 1L,
                                    malignancy=='no' ~ 0L)) %>%
    group_by(person_id) %>%
    summarise(progressive_ct = sum(as.numeric(progressive_num)),
              malignancy_ct=sum(as.numeric(malignancy_num)))
  
  body_systems_count <- pmca_summary_tbl %>%
    distinct(person_id, body_system) %>%
    mutate(body_system_flag=1L) %>%
    ungroup() %>%
    group_by(person_id) %>%
    summarise(n_body_systems = count(body_system_flag)) %>%
    ungroup()
  
  
  final_flags <- pmca_summary_tbl %>%
    distinct(person_id) %>%
    left_join(body_systems_count) %>%
    left_join(progressive_malignancy_flag) %>%
    right_join(select(cohort_tbl, person_id)) %>%
    right_join(select(pmca_category, person_id, site, complex_chronic, chronic, non_complex_chronic)) %>%
    mutate(
      progressive_ct = case_when(progressive_ct > 0L ~ progressive_ct,
                                 TRUE ~ 0L),
      malignancy_ct = case_when(malignancy_ct > 0L ~ malignancy_ct,
                                TRUE ~ 0L),
      n_body_systems = case_when(n_body_systems > 0L ~ n_body_systems,
                                 TRUE ~ 0L)) %>%
    select(person_id,
           site,
           progressive_ct,
           malignancy_ct,
           complex_chronic,
           chronic,
           non_complex_chronic,
           n_body_systems)
  
}





#' Function from pasc25_hazard_ratio / code / cohort_pmca -- modified to produce by site
#' produce PMCA table with 3 year lookback
#'
#' @param cohort cohort of patients with patients and observation_date
#' @param pmca_xwalk codeset that has flags for body systems and whether or not progressive
#' @param condition_tbl formatting of the condition occurrence table
#' @param visit_tbl formatting of the visit occurrence table
#'
#' @return table that has conditions and flags for body systems
#' columns:
#' person_id | observation_date | result_derivation | condition_concept_id | condition_concept_name |
#' condition_source_concept_id | condition_source_value | visit_occurrence_id | condition_start_state |
#' description | body_system | progressive | visit_concept_id | site
#'

produce_pmca_lookup_by_site <- function (cohort,
                                         pmca_xwalk=load_codeset('pmca_icd10'),
                                         condition_tbl=cdm_tbl('condition_occurrence'),
                                         visit_tbl=cdm_tbl('visit_occurrence'),
                                         site='chop') {
  
  cond_tbl <- condition_tbl %>%
    filter(site==site)
  
  visit_occ_tbl <- visit_tbl %>%
    filter(site==site)
  
  only_pmca_conds <-
    cohort %>%
    inner_join(
      select(cond_tbl,
             person_id,
             condition_concept_id,
             condition_concept_name,
             condition_source_concept_id,
             condition_source_concept_name,
             visit_occurrence_id,
             condition_start_date),
      by=c('person_id')
    ) %>% filter(observation_date - condition_start_date <= 1096L &
                   observation_date - condition_start_date >= 0L) %>%
    
    inner_join(
      pmca_xwalk,
      by=c('condition_source_concept_id'='concept_id')
    ) %>%
    inner_join(
      select(
        visit_occ_tbl,
        visit_occurrence_id,
        visit_concept_id
      )
    ) %>% filter(visit_concept_id %in% c(9202L,9201L,9203L,
                                         2000000048L,2000000088L,581399L)) %>%
    distinct() 
  
}




compute_pmca_site <- function(cohort,
                              site_list=list('chop','cchmc',
                                             'colorado','lurie',
                                             'nationwide','nemours',
                                             'seattle','stanford')) {
  
  site_rslts <- list()
  
  for(i in 1:length(site_list)) {
    
    site_nm = site_list[[i]]
    
    pmca_lookup <- produce_pmca_lookup_by_site(cohort=cohort %>% filter(site==site_nm),
                                               site=site_nm)
    
    pmca_summary <- compute_pmca_summary(pmca_lookup_tbl = pmca_lookup) %>%
      mutate(site=site_nm)
    
    site_rslts[[i]] <- pmca_summary
    
  }
  
  site_rslts
  
}





#' Flag PMCA body systems, per system
#' @param cohort_tbl cohort table with list of patients
#' @param pmca_summary pmca summary tbl
#' @return a table with body systems listed
#' 
flag_pmca_bs <- function(cohort_tbl=results_tbl('cohort_all', db=config('db_src_full')),
                         pmca_summary=results_tbl('pmca_summary_full', db=config('db_src_full'))) {
  
  pmca_bs <- pmca_summary %>%
    filter(total_visits>0L) %>%
    distinct(person_id, body_system) %>%
    mutate(has_body_system=1L) %>%
    pivot_wider(names_from='body_system', values_from='has_body_system', values_fill=0) %>%
    full_join(cohort_tbl %>% distinct(person_id), by='person_id') %>%
    rename(pulmonary_respiratory=`pulmonary/respiratory`,
           mental_health=`mental health`) %>%
    mutate(hematological = case_when(is.na(hematological) ~ 0L,
                                     TRUE ~ hematological),
           renal = case_when(is.na(renal) ~ 0L,
                             TRUE ~ renal),
           gastrointestinal = case_when(is.na(gastrointestinal) ~ 0L,
                                        TRUE ~ gastrointestinal),
           musculoskeletal = case_when(is.na(musculoskeletal) ~ 0L,
                                       TRUE ~ musculoskeletal),
           pulmonary_respiratory = case_when(is.na(pulmonary_respiratory) ~ 0L,
                                             TRUE ~ pulmonary_respiratory),
           genitourinary = case_when(is.na(genitourinary) ~ 0L,
                                     TRUE ~ genitourinary),
           otologic = case_when(is.na(otologic) ~ 0L,
                                TRUE ~ otologic),
           genetic = case_when(is.na(genetic) ~ 0L,
                               TRUE ~ genetic),
           metabolic = case_when(is.na(metabolic) ~ 0L,
                                 TRUE ~ metabolic),
           endocrinological = case_when(is.na(endocrinological) ~ 0L,
                                        TRUE ~ endocrinological),
           craniofacial = case_when(is.na(craniofacial) ~ 0L,
                                    TRUE ~ craniofacial),
           cardiac = case_when(is.na(cardiac) ~ 0L,
                               TRUE ~ cardiac),
           immunological = case_when(is.na(immunological) ~ 0L,
                                     TRUE ~ immunological),
           dermatological = case_when(is.na(dermatological) ~ 0L,
                                      TRUE ~ dermatological),
           ophthalmological = case_when(is.na(ophthalmological) ~ 0L,
                                        TRUE ~ ophthalmological),
           mental_health = case_when(is.na(mental_health) ~ 0L,
                                     TRUE ~ mental_health),
           neurological = case_when(is.na(neurological) ~ 0L,
                                    TRUE ~ neurological)) %>%
    compute_new(indexes=c('person_id'))
  
}








