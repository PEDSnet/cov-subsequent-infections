
.run <- function(){
  
  # expand codeset for high_risk conditions
  # dx_highrisk <- load_codeset("dx_highrisk") %>% collect()
  # highrisk_codeset <- list()
  # for (i in 1:dim(dx_highrisk)[1]){
  #   temp_codes <- dx_highrisk[i, "concept_code"] %>% pull()
  #   temp <- vocabulary_tbl("concept") %>%
  #     filter(grepl(temp_codes, concept_code, ignore.case = TRUE),
  #            grepl("ICD10", vocabulary_id, ignore.case = TRUE)) %>%
  #     select(concept_id, concept_name, concept_code, vocabulary_id) %>% collect()
  #   highrisk_codeset <- rbind.data.frame(highrisk_codeset, temp)
  # }
  # highrisk_codeset %>% view()
  # highrisk_codeset %>% output_tbl("dx_highrisk_final", file = TRUE, local = TRUE)
  
  # get analytic dataset:
  # postacute_min_start_date = as.Date("2022-04-01")
  # postacute_max_end_date = as.Date("2023-01-01")
  # 
  # cohort_entry_start_date = as.Date("2022-03-01")
  # cohort_entry_end_date = as.Date("2022-07-01")
  # 
  # cohort_1_label = "rsv_study_cohort"
  
  postacute_min_start_date_resp = as.Date("2021-07-01")
  postacute_max_end_date_resp = as.Date("2023-03-01")
  
  cohort_entry_start_date_resp = as.Date("2021-06-01")
  cohort_entry_end_date_resp = as.Date("2023-02-01")
  
  cohort_1_label = "resp_study_cohort"
  comparison_group_string = "Influenza"
  
  comparison_group_string = "Influenza"
  
  dx_highrisk <- load_codeset("dx_highrisk")
  
  analytic_dataset <-
    results_tbl(paste0(cohort_1_label, "_analytic_dataset")) %>% compute_new()
  # highrisk_tbl <- analytic_dataset %>% slice_sample(n = 10) %>% select(person_id)
  
  highrisk_tbl <- analytic_dataset %>% select(person_id) %>% 
    left_join(cdm_tbl("condition_occurrence") %>% 
                select(condition_source_concept_id, person_id),# %>% compute_new(indexes=c("person_id", "condition_source_concept_id")), 
              by="person_id") %>% 
    distinct(person_id, condition_source_concept_id, .keep_all = TRUE)%>%
    compute_new(indexes=c("person_id", "condition_source_concept_id"))
  highrisk_tbl <- highrisk_tbl %>% inner_join(dx_highrisk %>% select(concept_id) %>% mutate(highrisk_flag = TRUE),
                                              by = c("condition_source_concept_id" = "concept_id")) %>%
    select(-condition_source_concept_id) %>%
    distinct()%>% compute_new()
  
  analytic_dataset <- analytic_dataset %>% left_join(highrisk_tbl, by = "person_id") %>%
    mutate(highrisk_flag = if_else(is.na(highrisk_flag), FALSE, TRUE))
  analytic_dataset %>% output_tbl(paste0(cohort_1_label, "_analytic_dataset_highrisk"))
}
