## Codesets for flu

### Flu diagnosis

flu_condition_codes <-
  vocabulary_tbl('concept') %>% 
  filter(str_detect(concept_code, 'J09') | 
           str_detect(concept_code, 'J10') | 
           str_detect(concept_code, 'J11')) %>% 
  filter(vocabulary_id %in% c("ICD10", "ICD10CM"))

flu_conditions_snomed <-
  vocabulary_tbl("concept_relationship") %>% 
  inner_join(flu_condition_codes, by=c("concept_id_1"="concept_id")) %>% 
  filter(relationship_id=="Maps to") %>% 
  select(concept_id = concept_id_2) %>% 
  left_join(vocabulary_tbl("concept"), by=c("concept_id"))

flu_condition_codes %>% 
  select(concept_id, concept_name, concept_code, vocabulary_id) %>% 
  full_join(flu_conditions_snomed %>% 
              select(concept_id, concept_name, concept_code, vocabulary_id)) %>% 
  write.csv('specs/dx_flu.csv', row.names = FALSE)



## Codeset for bacterial respiratory infections
# BACTERIA: 
#   
# ● Bordetella parapertussis 
# ● Bordetella pertussis 
# ● Chlamydia pneumoniae 
# ● Mycoplasma pneumoniae 

bacterial_resp_concepts <- vocabulary_tbl('concept') %>% 
  filter(str_detect(tolower(concept_name), 'bordatella parapertussis') |
           str_detect(tolower(concept_name), 'mycoplasma pneumoniae') |
           str_detect(tolower(concept_name), 'chlamydia pneumoniae') |
           str_detect(tolower(concept_name), 'bordetella pertussis')) %>% 
  filter(vocabulary_id %in% c("SNOMED", "ICD10", "ICD10-CM")) %>% 
  filter(domain_id %in% c('Condition')) %>% 
  compute_new(indexes=c("concept_id"))

bacterial_resp_concepts %>% 
  select(concept_id, concept_name, concept_code, vocabulary_id) %>% 
  write.csv('specs/bacterial_respiratory_infections_draft.csv', row.names = FALSE)

### Influenza ILI
## TODO: Get descendants
## TODO: Map to ICD10 using concept relationship
# influenza_and_ili_conditions <- vocabulary_tbl('concept') %>% 
#   filter(str_detect(tolower(concept_name), 'influenza')) %>% 
#   filter(!(str_detect(tolower(concept_name), 'vaccination'))) %>% 
#   filter(!(str_detect(tolower(concept_name), 'vaccine'))) %>% 
#   filter(vocabulary_id %in% c("SNOMED", "ICD10", "ICD10-CM")) %>% 
#   filter(domain_id %in% c('Condition')) %>% 
#   compute_new(indexes=c("concept_id"))
# 
# influenza_and_ili_conditions %>% 
#   select(concept_id, concept_name, concept_code, vocabulary_id) %>% 
#   write.csv('specs/ili_conditions.csv', row.names = FALSE)

## Any infection
## Infectious disease, aka "disorder due to infection" (synonym)

infectious_disease_snomed <- vocabulary_tbl('concept') %>% 
  filter(vocabulary_id %in% c("SNOMED-CT", "SNOMED"))  %>% 
  filter(concept_code == "40733004") %>% 
  compute_new(indexes=c("concept_id")) %>% 
  get_descendants() %>% 
  compute_new()

infectious_disease_snomed_filtered <-
  infectious_disease_snomed %>% 
  filter(vocabulary_id %in% c("SNOMED-CT", "SNOMED"))

snomed_to_icd10_infection_codes <-
  vocabulary_tbl("concept_relationship") %>% 
  inner_join(infectious_disease_snomed_filtered, by=c("concept_id_1"="concept_id")) %>% 
  filter(relationship_id=="Mapped from") %>% 
  select(concept_id = concept_id_2) %>% 
  left_join(vocabulary_tbl("concept"), by=c("concept_id")) %>% 
  filter(vocabulary_id %in% c("ICD10", "ICD10CM")) %>% 
  compute_new()


infectious_disease_snomed_filtered %>% 
  select(concept_id, concept_name, concept_code, vocabulary_id) %>% 
  full_join(snomed_to_icd10_infection_codes %>% 
              select(concept_id, concept_name, concept_code, vocabulary_id), 
            by=c("concept_id", "concept_name", "concept_code", "vocabulary_id")) %>% 
  write.csv('specs/dx_any_infection.csv', row.names = FALSE)

### Excluding covid from any infection codeset:
any_infection_codeset_without_covid <- load_codeset("dx_any_infection") %>% 
  filter(!(
    str_detect(tolower(concept_name), 'covid') |
      str_detect(tolower(concept_name), 'sars-cov') |
      str_detect(tolower(concept_name), 'coronavirus')))

any_infection_codeset_without_covid %>% 
  write.csv('specs/dx_noncov_any_infection.csv', row.names = FALSE)



### Lab_rsv

lab_rsv_strings <- vocabulary_tbl('concept') %>% 
  filter(str_detect(tolower(concept_name), 'rsv') |
           str_detect(tolower(concept_name), 'respiratory syncytial virus')) %>% 
  filter(vocabulary_id %in% c("LOINC")) %>% 
  compute_new(indexes=c("concept_id")) %>% 
  get_descendants() %>% 
  compute_new()

lab_rsv_strings %>% 
  select(concept_id, concept_name, concept_code, vocabulary_id) %>% 
  write.csv('specs/lab_rsv.csv', row.names = FALSE)


### Dx RSV, subset of dx_respiratory infections

respiratory_codeset <- load_codeset("dx_noncov_resp_infections")

rsv_dx <- respiratory_codeset %>% 
  filter(str_detect(tolower(concept_name), 'rsv') |
           str_detect(tolower(concept_name), 'respiratory syncytial virus'))

## TODO: needs to be filtered for whether or not the "due to" category is applicable
rsv_dx %>% 
  select(concept_id, concept_name, concept_code, vocabulary_id) %>% 
  write.csv('specs/dx_rsv_draft.csv', row.names = FALSE)

# ## Lab MSK
# lab_rsv_strings <- vocabulary_tbl('concept') %>% 
#   filter(str_detect(tolower(concept_name), 'rsv') |
#            str_detect(tolower(concept_name), 'respiratory syncytial virus')) %>% 
#   filter(vocabulary_id %in% c("LOINC")) %>% 
#   compute_new(indexes=c("concept_id")) %>% 
#   get_descendants() %>% 
#   compute_new()
# 
# lab_rsv_strings %>% 
#   select(concept_id, concept_name, concept_code, vocabulary_id) %>% 
#   write.csv('specs/lab_rsv.csv', row.names = FALSE)
# 
# ## Lab SSTI
# ssti_strings <- vocabulary_tbl('concept') %>% 
#   filter(str_detect(tolower(concept_name), 'skin') |
#            str_detect(tolower(concept_name), 'soft tissue')) %>% 
#   filter(vocabulary_id %in% c("LOINC")) %>% 
#   compute_new(indexes=c("concept_id")) %>% 
#   get_descendants() %>% 
#   compute_new()


#### Modifying existing codesets
rsv_codeset_with_sr_edits <- load_codeset("lab_rsv")

filtered_lab_rsv <- rsv_codeset_with_sr_edits %>% 
  filter(is.na(exclude) | exclude != "x") 

filtered_lab_rsv %>% 
  select(concept_id, concept_name, concept_code, vocabulary_id) %>% 
  write.csv('specs/lab_rsv_final.csv', row.names = FALSE)

## Updates: don't use dx_influenza, now just use resp codes filtered down.
## Use lab_rsv_final for RSV outcome codes
## use resp codes filtered down to the RSV codes for rsv

## now we have covid, influenza, and all other resp cohort (non-influenza and non-rsv)
## and then we have rsv dx's and labs as outcomes

load_codeset("dx_noncov_resp_infections") %>% 
  filter(influenza_related == "x") %>% 
  select(concept_id, concept_name, concept_code, vocabulary_id) %>% 
  write.csv('specs/dx_influenza_final.csv', row.names = FALSE)

load_codeset("dx_noncov_resp_infections") %>% 
  filter(rsv_related == "x") %>% 
  select(concept_id, concept_name, concept_code, vocabulary_id) %>% 
  write.csv('specs/dx_rsv_final.csv', row.names = FALSE)

load_codeset("dx_noncov_resp_infections") %>% 
  filter(is.na(rsv_related), is.na(influenza_related)) %>% 
  select(concept_id, concept_name, concept_code, vocabulary_id) %>% 
  write.csv('specs/dx_other_resp_final.csv', row.names = FALSE)

### FINAL CODESETS
## dx_other_resp_final
## dx_rsv_final
## dx_influenza_final
## lab_rsv_final
## lab_influenza_filtered


  
