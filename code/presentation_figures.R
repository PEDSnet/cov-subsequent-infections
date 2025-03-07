exclusion_reasons <- results_tbl(paste0(cohort_1_label, "_w_exclusion_reasons"))

attrition <- results_tbl(paste0(cohort_1_label, "_attrition"))

final_cohort_demo <- results_tbl(paste0(cohort_1_label, "_cohort_demo"))

rsv_outcomes <- results_tbl(paste0(cohort_1_label, "_rsv_outcomes"))


cohort_with_rsv_outcomes <-
  final_cohort_demo %>% 
  left_join(rsv_outcomes, by="person_id")


rsv_week_evidence <- cohort_with_rsv_outcomes %>% 
  mutate(rsv_evidence_week = floor_date(rsv_evidence_date, unit="week")) %>% 
  mutate(rsv_evidence_week = as.Date(rsv_evidence_week)) %>% 
  group_by(rsv_evidence_week) %>% 
  summarise(n_pats = n_distinct(person_id)) %>% 
  select(rsv_evidence_week, n_pats) %>% 
  filter(!is.na(rsv_evidence_week)) %>% 
  compute_new()

resp_week_evidence <- cohort_with_rsv_outcomes %>% 
  mutate(resp_week = floor_date(resp_date, unit="week")) %>% 
  mutate(resp_evidence_week = as.Date(resp_week)) %>% 
  group_by(resp_evidence_week) %>% 
  summarise(n_pats = n_distinct(person_id)) %>% 
  select(resp_evidence_week, n_pats) %>% 
  filter(!is.na(resp_evidence_week)) %>% 
  compute_new()

flu_week_evidence <- cohort_with_rsv_outcomes %>% 
  mutate(flu_week = floor_date(flu_date, unit="week")) %>% 
  mutate(flu_week = as.Date(flu_week)) %>% 
  group_by(flu_week) %>% 
  summarise(n_pats = n_distinct(person_id)) %>% 
  select(flu_week, n_pats) %>% 
  filter(!is.na(flu_week)) %>% 
  compute_new()

covid_week_evidence <- cohort_with_rsv_outcomes %>% 
  mutate(covid_week = floor_date(covid_date, unit="week")) %>% 
  mutate(covid_evidence_week = as.Date(covid_week)) %>% 
  group_by(covid_evidence_week) %>% 
  summarise(n_pats = n_distinct(person_id)) %>% 
  select(covid_evidence_week, n_pats) %>% 
  filter(!is.na(covid_evidence_week)) %>% 
  compute_new()

all_weekly_incidence_data <-
  covid_week_evidence %>% 
  select(week = covid_evidence_week, covid_pats = n_pats) %>% 
  full_join(flu_week_evidence %>% 
              select(week = flu_week, flu_pats = n_pats), by="week") %>% 
  full_join(rsv_week_evidence %>% 
              select(week = rsv_evidence_week, rsv_pats = n_pats), by="week") %>% 
  full_join(resp_week_evidence %>% 
              select(week = resp_evidence_week, resp_pats = n_pats), by="week") %>% 
  compute_new()

all_weekly_incidence_data %>% 
  pivot_longer(cols=c("covid_pats", "flu_pats", "rsv_pats", "resp_pats"),
               names_to = "disease",
               values_to = "n_pats") %>% 
  filter(!is.na(n_pats)) %>% 
  ggplot() +
  geom_point(aes(x=week, y = n_pats, color = disease), alpha=0.5) +
  geom_smooth(aes(x=week, y = n_pats, color = disease), alpha=0.5, span = .2)

cohort_with_rsv_outcomes %>% 
  ggplot() +
  geom_density(aes(x=ce_date, fill = sub_cohort), alpha=0.2, position="stack") +
  geom_density(aes(x=rsv_evidence_date, fill = sub_cohort), alpha=0.2, position="stack") +
  theme_bw()

## Exploring the outcome of RSV vs. the index events

cohort_with_rsv_outcomes %>% 
  filter(!is.na(rsv_evidence_date)) %>% 
  ggplot() +
  geom_density(aes(x=ce_date, fill = sub_cohort), alpha=0.2, position="stack") +
  geom_density(aes(x=rsv_evidence_date, fill = "rsv date"), alpha=0.2, position="stack") +
  theme_bw() +
  facet_wrap(~sub_cohort)

cohort_with_rsv_outcomes %>% 
  filter(!is.na(rsv_evidence_date)) %>% 
  mutate(days_from_ce_to_rsv = as.numeric(rsv_evidence_date - ce_date)) %>% 
  ggplot() +
  geom_density(aes(x=days_from_ce_to_rsv, fill = sub_cohort), alpha=0.2, position="stack") +
  theme_bw() +
  facet_wrap(~sub_cohort)





