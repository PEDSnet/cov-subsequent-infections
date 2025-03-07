output_attrition_table <- function(cohort_1_label, cohort_entry_start_date, cohort_entry_end_date){
  init_sum(cohort = 'Start', persons = 0)
  
  append_sum(cohort = 'Patients under 6 years old with at least 1 visit during CE period',
             persons = distinct_ct(results_tbl(paste0(cohort_1_label, "_pats_under_6"))))
  if(cohort_1_label == "resp_study"){
    append_sum(cohort = 'Had at least 1 of COVID, influenza infection event during CE period',
               persons = distinct_ct(results_tbl(paste0(cohort_1_label, "_triple_cohort"))))
  } else {
    append_sum(cohort = 'Had at least 1 of COVID, influenza, and/or other respiratory infection event during CE period',
               persons = distinct_ct(results_tbl(paste0(cohort_1_label, "_triple_cohort"))))
  }
  
  append_sum(cohort = 'Had COVID during CE period',
             persons = distinct_ct(results_tbl(paste0(cohort_1_label, "_triple_cohort")) %>% 
                                     filter(!is.na(covid_date),
                                            covid_date >= cohort_entry_start_date,
                                            covid_date < cohort_entry_end_date)))
  append_sum(cohort = 'Had influenza during CE period',
             persons = distinct_ct(results_tbl(paste0(cohort_1_label, "_triple_cohort")) %>%
                                     filter(!is.na(flu_date),
                                            flu_date >= cohort_entry_start_date,
                                            flu_date < cohort_entry_end_date)))
  if(cohort_1_label == "rsv_study"){
    append_sum(cohort = 'Had other respiratory infection during CE period',
               persons = distinct_ct(results_tbl(paste0(cohort_1_label, "_triple_cohort")) %>%
                                       filter(!is.na(resp_date),
                                              resp_date >= cohort_entry_start_date,
                                              resp_date < cohort_entry_end_date)))
  }
  ## TODO add to attrition step (study_eligible == 1)
  append_sum(cohort = 'Meets utilization inclusion criteria',
             persons = distinct_ct(results_tbl(paste0(cohort_1_label, "_cohort_inclusion_flags")) %>% 
                                     filter(study_eligible == 1)))
  if(cohort_1_label == "resp_study"){
    append_sum(cohort = 'Overall cohort without overlapping index infections',
               persons = distinct_ct(results_tbl(paste0(cohort_1_label, "_w_excl_reasons")) %>% 
                                       filter(study_eligible == 1,
                                              exclude != 1)))
  } else {
    append_sum(cohort = 'Overall cohort without overlapping index infections',
               persons = distinct_ct(results_tbl(paste0(cohort_1_label, "_other_washout_reasons")) %>% 
                                       filter(study_eligible == 1,
                                              exclude != 1)))
  }
  
  
  rslt$final_cohort <- results_tbl(paste0(cohort_1_label, "_cohort"))
  append_sum(cohort = 'COVID sub-cohort',
             persons = distinct_ct(rslt$final_cohort %>% 
                                     filter(sub_cohort=="COVID")))
  
  append_sum(cohort = 'Influenza sub-cohort',
             persons = distinct_ct(rslt$final_cohort %>% 
                                     filter(sub_cohort=="Influenza")))
  if(cohort_1_label == "rsv_study"){
    append_sum(cohort = 'Other respiratory infection sub-cohort',
               persons = distinct_ct(rslt$final_cohort %>% 
                                       filter(sub_cohort=="Respiratory")))
  }
  
  analytic_dataset_final<- results_tbl(paste0(cohort_1_label, "_analytic_final"))%>%
    filter(ce_date >= cohort_entry_start_date,
           ce_date < as.Date("2023-01-01")) %>%
    filter(sex_cat != "Other/unknown") %>% 
    filter(is.na(covid_index_date_imputed) | covid_index_date_imputed==0)
  
  if(cohort_1_label == "rsv_study") {
    analytic_dataset_final <- analytic_dataset_final%>%
      group_by(person_id) %>% 
      mutate(exclude = max(exclude_for_prior_rsv)) %>% 
      filter(exclude < 1) %>%
      slice_min(rsv_evidence_date, with_ties=FALSE) %>% 
      ungroup() 
  }
  
  append_sum(cohort = 'Overall cohort without invalid gender and co-infections',
             persons = analytic_dataset_final %>% distinct_ct())
  append_sum(cohort = 'Final COVID sub-cohort',
             persons = analytic_dataset_final %>%
               filter(sub_cohort=="COVID") %>% distinct_ct())
  append_sum(cohort = 'Final Influenza sub-cohort',
             persons = analytic_dataset_final %>%
               filter(sub_cohort=="Influenza") %>% distinct_ct())
  if(cohort_1_label == "rsv_study") {
    append_sum(cohort = 'Final Respiratory sub-cohort',
               persons = analytic_dataset_final %>%
                 filter(sub_cohort=="Respiratory") %>% distinct_ct())
  }
  output_sum(name = paste0(cohort_1_label, "_attrition"))
  
  return()
}

# exclude certain sites
output_attrition_table_site_excl <- function(cohort_1_label, cohort_entry_start_date, cohort_entry_end_date, site_excl){
  
  init_sum(cohort = 'Start', persons = 0)
  
  append_sum(cohort = 'Patients under 6 years old with at least 1 visit during CE period',
             persons = distinct_ct(results_tbl(paste0(cohort_1_label, "_pats_under_6")) %>% exclude_sites(site_excl)))
  if(cohort_1_label == "resp_study"){
    append_sum(cohort = 'Had at least 1 of COVID, influenza infection event during CE period',
               persons = distinct_ct(results_tbl(paste0(cohort_1_label, "_triple_cohort")) %>% exclude_sites(site_excl)))
  } else {
    append_sum(cohort = 'Had at least 1 of COVID, influenza, and/or other respiratory infection event during CE period',
               persons = distinct_ct(results_tbl(paste0(cohort_1_label, "_triple_cohort")) %>% exclude_sites(site_excl)))
  }
  
  append_sum(cohort = 'Had COVID during CE period',
             persons = distinct_ct(results_tbl(paste0(cohort_1_label, "_triple_cohort")) %>% 
                                     filter(!is.na(covid_date),
                                            covid_date >= cohort_entry_start_date,
                                            covid_date < cohort_entry_end_date) %>% exclude_sites(site_excl)))
  append_sum(cohort = 'Had influenza during CE period',
             persons = distinct_ct(results_tbl(paste0(cohort_1_label, "_triple_cohort")) %>%
                                     filter(!is.na(flu_date),
                                            flu_date >= cohort_entry_start_date,
                                            flu_date < cohort_entry_end_date) %>% exclude_sites(site_excl)))
  if(cohort_1_label == "rsv_study"){
    append_sum(cohort = 'Had other respiratory infection during CE period',
               persons = distinct_ct(results_tbl(paste0(cohort_1_label, "_triple_cohort")) %>%
                                       filter(!is.na(resp_date),
                                              resp_date >= cohort_entry_start_date,
                                              resp_date < cohort_entry_end_date) %>% exclude_sites(site_excl)))
  }
  ## TODO add to attrition step (study_eligible == 1)
  append_sum(cohort = 'Meets utilization inclusion criteria',
             persons = distinct_ct(results_tbl(paste0(cohort_1_label, "_cohort_inclusion_flags")) %>% 
                                     filter(study_eligible == 1) %>% exclude_sites(site_excl)))
  if(cohort_1_label == "resp_study"){
    append_sum(cohort = 'Overall cohort without overlapping index infections',
               persons = distinct_ct(results_tbl(paste0(cohort_1_label, "_w_excl_reasons")) %>% 
                                       filter(study_eligible == 1,
                                              exclude != 1) %>% exclude_sites(site_excl)))
  } else {
    append_sum(cohort = 'Overall cohort without overlapping index infections',
               persons = distinct_ct(results_tbl(paste0(cohort_1_label, "_other_washout_reasons")) %>% 
                                       filter(study_eligible == 1,
                                              exclude != 1) %>% exclude_sites(site_excl)))
  }
  
  
  rslt$final_cohort <- results_tbl(paste0(cohort_1_label, "_cohort"))  %>% exclude_sites(site_excl)
  append_sum(cohort = 'COVID sub-cohort',
             persons = distinct_ct(rslt$final_cohort %>% 
                                     filter(sub_cohort=="COVID") ))
  
  append_sum(cohort = 'Influenza sub-cohort',
             persons = distinct_ct(rslt$final_cohort %>% 
                                     filter(sub_cohort=="Influenza")) )
  if(cohort_1_label == "rsv_study"){
               append_sum(cohort = 'Other respiratory infection sub-cohort',
                          persons = distinct_ct(rslt$final_cohort %>% 
                                                  filter(sub_cohort=="Respiratory") ))
             }
             
             analytic_dataset_final<- results_tbl(paste0(cohort_1_label, "_analytic_final"))%>%
               filter(ce_date >= cohort_entry_start_date,
                      ce_date < as.Date("2023-01-01")) %>%
               filter(sex_cat != "Other/unknown") %>% 
               filter(is.na(covid_index_date_imputed) | covid_index_date_imputed==0) %>% 
               exclude_sites(site_excl)
             
             if(cohort_1_label == "rsv_study") {
               analytic_dataset_final <- analytic_dataset_final%>%
                 group_by(person_id) %>% 
                 mutate(exclude = max(exclude_for_prior_rsv)) %>% 
                 filter(exclude < 1) %>%
                 slice_min(rsv_evidence_date, with_ties=FALSE) %>% 
                 ungroup() 
             }
             
             append_sum(cohort = 'Overall cohort without invalid gender and co-infections',
                        persons = analytic_dataset_final %>% distinct_ct())
             append_sum(cohort = 'Final COVID sub-cohort',
                        persons = analytic_dataset_final %>%
                          filter(sub_cohort=="COVID") %>% distinct_ct())
             append_sum(cohort = 'Final Influenza sub-cohort',
                        persons = analytic_dataset_final %>%
                          filter(sub_cohort=="Influenza") %>% distinct_ct())
             if(cohort_1_label == "rsv_study") {
               append_sum(cohort = 'Final Respiratory sub-cohort',
                          persons = analytic_dataset_final %>%
                            filter(sub_cohort=="Respiratory") %>% distinct_ct())
             }
             output_sum(name = paste0(cohort_1_label, "_attrition_site_excl"))
             
             return()
}

exclude_sites <- function(cohort, site_exclude_list){
  cohort %>% inner_join(cdm_tbl("person") %>% select(person_id, site), by = "person_id") %>%
    filter(!site %in% site_exclude_list) %>%
    return()
}
