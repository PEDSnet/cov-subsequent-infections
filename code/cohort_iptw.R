
library(ggpubr)

do_iptw <- function(analytic_dataset, cohort_entry_start_date, comparison_cohort,
                    weight_formula = 'covid~age_group+sex_cat+race_eth_cat+ce_week_factor', 
                    iptw_description) {
  #### IPTW
  
  iptw_dataset <-
    analytic_dataset %>% 
    mutate(ce_week = as.Date(floor_date(ce_date, unit="week"))) %>% 
    mutate(study_start_date = as.Date(cohort_entry_start_date)) %>% 
    mutate(ce_date_days_numeric = as.numeric(ce_date-study_start_date)) %>% 
    mutate(ce_week_days_numeric = as.numeric(ce_week-study_start_date)) %>% 
    collect() %>% 
    mutate(ce_week_factor = as.factor(ce_week_days_numeric)) 
  # compute_new(indexes=c("person_id"))
  
  ### Now first want to limit the cohort to the COVID/flu comparison:
  iptw_comparison_groups <-
    iptw_dataset %>% 
    filter(sub_cohort %in% c("COVID", comparison_cohort)) %>% 
    mutate(covid = ifelse(sub_cohort=="COVID", 1, 0)) 
  
  
  #### Calling modeling stuff from function
  # weight_formula = as.formula('covid~age_group+sex_cat+race_eth_cat+ce_week_factor')
  
  # weight_formula_with_util = as.formula('covid~age_group+sex_cat+race_eth_cat+ce_week_factor+util_other+util_inpatient+util_outpatient+util_ed')
  
  # TODO change the function also to output model coefficient results based on pre-weight analysis
  iptw_results <- 
    make_iptw_and_mw_weights_colby(weight_formula = as.formula(weight_formula),
                                   dataset = iptw_comparison_groups,
                                   description = iptw_description # covid_flu_noutil
    )
  
  iptw_results %>% 
    return()
}

#' Code from Katie Hirabayashi

#' Use WeightIt package to make weights for IPTW; and add MW
#' @param weight_formula formula to weight your model; make sure outcome is entered as binary 1/0
#' @param dataset dataset for weighting
#' @param description string descriptor for output
#' @return dataset with propensity scores, IPTW weights, IPTW & matching weights (MW), and 95th percentile-capped IPTW weights
#' samplesize (sample size), weightrange (highest and lowest weights), weighttop (highest weights),
#' coefofvar (coefficient of variation), iptw_blance (balancetable for iptw), iptw_observations (iptw obserations),
#' iptw (dataset with capped weights?) maxweights_postcap (maximum weights post 95% weight capping)
make_iptw_and_mw_weights_colby <- function(weight_formula=as.formula('covid~age_group+
                                                                       sex_cat+
                                                                       raceth_cat+
                                                                       total_pre_1y+
                                                                       complex_chronic_flag+
                                                                       prior_sud+
                                                                       bmi_cat+
                                                                       glucocorticoids+
                                                                       plan_class+
                                                                       index_period'),
                                           dataset=rslt$cov_vs_oriac,
                                           description='iptw_covid_oriac') {
  
  dataset_final <- set_iptw_levels(iptw_data=dataset)
  
  #' Propensity Score Model
  ps_model <- glm(weight_formula,
                  family=binomial(link="logit"),
                  data=dataset_final)
  
  pscore <- ps_model$fitted.values
  
  dataset_final$p_score <- predict(ps_model, type = "response")
  
  visuals = list()
  viz_pscores <- dataset_final %>% 
    select(person_id, sub_cohort, covid, p_score) %>% 
    ggplot() + 
    geom_boxplot(aes(x=sub_cohort, y=p_score)) +
    theme_bw() +
    scale_color_brewer(palette = "Dark2") +
    theme(legend.position = "top")
  
  visuals$uncapped_pscores <- viz_pscores
  
  lower_cut <- dataset_final %>% filter(covid==1L) %>% pull(p_score) %>% quantile(0.005)
  upper_cut <- dataset_final %>% filter(covid==0L) %>% pull(p_score) %>% quantile(0.995)
  
  dataset_final_trim <- dataset_final %>%
    filter((p_score > lower_cut) & (p_score < upper_cut))
  
  ps_model_trim <- glm(weight_formula,
                       family=binomial,
                       data=dataset_final_trim)
  
  pscore2 <- ps_model_trim$fitted.values
  
  dataset_final_trim$p_score <- predict(ps_model_trim, type = "response")
  
  viz_pscores_trim <- dataset_final_trim %>% 
    select(person_id, sub_cohort, covid, p_score) %>% 
    ggplot() + 
    geom_boxplot(aes(x=sub_cohort, y=p_score)) +
    theme_bw() +
    scale_color_brewer(palette = "Dark2") +
    theme(legend.position = "top")
  
  visuals$capped_pscores <- viz_pscores_trim
  
  
  #' Visualize overlap of raw propensity scores
  viz_ps_overlap <- ggplot(data=dataset_final_trim, aes(x=p_score, group=sub_cohort, fill=sub_cohort)) +
    geom_density(adjust=1.5, alpha=.4) +
    theme_bw() +
    scale_fill_brewer(palette = "Dark2") +
    theme(legend.position = "top")
  
  visuals$trimmed_ps_overlap <- viz_ps_overlap
  
  
  #' Insert weights to get traditional IPTW for each participant (ATE)
  dataset_final_trim$capped_iptw = ifelse((dataset_final_trim$covid==1L),
                                          (1/dataset_final_trim$p_score),
                                          (1/(1-dataset_final_trim$p_score))
  )
  
  ps_model_balanced <- glm(weight_formula, family=quasibinomial(link="logit"),
                           data = dataset_final_trim,
                           weights = capped_iptw)
  # Now see the covariates are not statistically significant anymore
  
  viz_pscore_v_iptw <-
    dataset_final_trim %>% 
    select(person_id, sub_cohort, covid, p_score, capped_iptw) %>% 
    ggplot() + 
    geom_point(aes(x=p_score, y = capped_iptw, color=sub_cohort)) +
    theme_bw() +
    scale_color_brewer(palette = "Dark2") +
    theme(legend.position = "top")
  
  visuals$trimmed_pscore_v_iptw <- viz_pscore_v_iptw
  
  
  #' Make a vector of the most extreme weights
  extreme_weights <- c('covid propensity score pre-trimming',
                       'control propensity score pre-trimming',
                       'covid propensity score post-trimming',
                       'control propensity score post-trimming',
                       'covid IPTW weight post-trimming',
                       'control IPTW weight post-trimming')
  
  extreme1 <- dataset_final %>% filter(covid==1L) %>% summarise(max_p=max(p_score)) %>% pull()
  extreme2 <- dataset_final %>% filter(covid==1L) %>% summarise(min_p=min(p_score)) %>% pull()
  extreme3 <- dataset_final %>% filter(covid==0L) %>% summarise(max_p=max(p_score)) %>% pull()
  extreme4 <- dataset_final %>% filter(covid==0L) %>% summarise(min_p=min(p_score)) %>% pull()
  extreme5 <- dataset_final_trim %>% filter(covid==1L) %>% summarise(max_p=max(p_score)) %>% pull()
  extreme6 <- dataset_final_trim %>% filter(covid==1L) %>% summarise(min_p=min(p_score)) %>% pull()
  extreme7 <- dataset_final_trim %>% filter(covid==0L) %>% summarise(max_p=max(p_score)) %>% pull()
  extreme8 <- dataset_final_trim %>% filter(covid==0L) %>% summarise(min_p=min(p_score)) %>% pull()
  extreme9 <- dataset_final_trim %>% filter(covid==1L) %>% summarise(max_iptw=max(capped_iptw)) %>% pull()
  extreme10 <- dataset_final_trim %>% filter(covid==1L) %>% summarise(min_iptw=min(capped_iptw)) %>% pull()
  extreme11 <- dataset_final_trim %>% filter(covid==0L) %>% summarise(max_iptw=max(capped_iptw)) %>% pull()
  extreme12 <- dataset_final_trim %>% filter(covid==0L) %>% summarise(min_iptw=min(capped_iptw)) %>% pull()
  extremes_high <- c(extreme1, extreme3, extreme5, extreme7, extreme9, extreme11)
  extremes_low <- c(extreme2, extreme4, extreme6, extreme8, extreme10, extreme12)
  
  extreme_final <- cbind(extreme_weights, extremes_high, extremes_low) %>%
    as_tibble() %>%
    mutate(extremes_high=as.numeric(extremes_high),
           extremes_low=as.numeric(extremes_low))
  
  #' Note that this function is producing the exact same weights as the manual insertion in the ifelse statement above
  #' (final_dataset_trimmed$capped_iptw and w=w.out1$weights are identical)
  #' But weightit has nice love plots and balance tables that can be used for visualization
  w.out1 <- WeightIt::weightit(weight_formula,
                               data = dataset_final_trim,
                               estimand = "ATE", method = "ps")
  
  set.cobalt.options(binary = "std")
  
  weights_auto <- summary(w.out1)
  
  final <- list()
  #final$ir_thou <- ir_thou

  final$love_plot <- love.plot(w.out1) +
    scale_color_brewer(palette = "Dark2") +
    theme(legend.position = "top")# gets the love plot weights using weightIt
  
  
  
  final$balance_tbl <- bal.tab(w.out1,
                               binary='std', continuous='std')
  final$smds <- final$balance_tbl$Balance %>% select(-Diff.Un) %>%
    rename(Diff_Adj=`Diff.Adj`) %>%
    mutate(SMD_gt_than_pt2=case_when(abs(Diff_Adj) > 0.2 ~ 1L,
                                     TRUE ~ 0L))
  final$extreme_final <- extreme_final
  
  ### GBM propensity score weighting
  # ps_model_gbm <- gbm(weight_formula,
  #                 distribution="bernoulli",
  #                 n.trees = 100,
  #                 interaction.depth = 4,
  #                 train.fraction = .8,
  #                 shrinkage = .0005,
  #                 data=dataset_final_trim)
  # 
  # visuals$viz_gbm <- summary(ps_model_gbm)
  # 
  # pscore_gbm <- ps_model_gbm$fitted.values
  # 
  # dataset_final_trim$p_score_gbm <- predict(ps_model_gbm, type = "response")
  # 
  # dataset_final_trim$iptw_gbm = ifelse((dataset_final_trim$covid==1L),
  #                                         (1/dataset_final_trim$p_score_gbm),
  #                                         (1/(1-dataset_final_trim$p_score_gbm))
  # )
  # 
  # viz_pscore_v_iptw_gbm <-
  #   dataset_final_trim %>% 
  #   select(person_id, sub_cohort, covid, p_score, p_score_gbm, capped_iptw, iptw_gbm) %>% 
  #   ggplot() + 
  #   geom_point(aes(x=p_score_gbm, y = iptw_gbm, color=sub_cohort), alpha=0.4) +
  #   theme_bw() +
  #   scale_color_brewer(palette = "Dark2") +
  #   theme(legend.position = "top") 
  # 
  # visuals$gbm_relative_influence <- visuals$viz_gbm %>% 
  #   ggplot() +
  #   geom_bar(aes(x=reorder(var, rel.inf), y = rel.inf, fill=rel.inf), stat="identity") +
  #   coord_flip() +
  #   theme_bw()

  
  final$dataset <- dataset_final_trim # gets the final dataset with the propensity scores using by-hand computation of the logistic regression
  
  
  # TODO find some way to save all the plots in one document or something
  visuals_plot_sheet <- ggarrange(visuals$uncapped_pscores, visuals$capped_pscores,
                                  visuals$trimmed_ps_overlap, visuals$trimmed_pscore_v_iptw,
                                  final$love_plot,
                                  # viz_pscore_v_iptw_gbm,
                                  # visuals$gbm_relative_influence,
                                  ncol = 2, nrow = 4, labels = "AUTO")
  ## Save the plot outcome as a summary sheet in results 
  ggsave(paste0("results/visual_summary_sheet_iptw_", description, ".png"), 
         plot = visuals_plot_sheet,
         width = 8.5, height = 11, dpi = 500
         )
  ggsave(paste0("results/visual_summary_sheet_iptw_", description, ".svg"), 
         plot = visuals_plot_sheet,
         width = 8.5, height = 11, dpi = 500
  )
  ggsave(paste0("results/visual_summary_sheet_iptw_", description, ".pdf"), 
         plot = visuals_plot_sheet,
         width = 8.5, height = 11, dpi = 500
  )
  
  #' Output final IPTW table with weights, post trimming
  output_tbl(dataset_final_trim, paste0(description,'_tabc'), indexes=c('person_id'))
  
  return(final)
  
}

library(gbm)
make_iptw_and_mw_weights_gbm <- function(weight_formula=as.formula('covid~age_group+
                                                                       sex_cat+
                                                                       raceth_cat+
                                                                       total_pre_1y+
                                                                       complex_chronic_flag+
                                                                       prior_sud+
                                                                       bmi_cat+
                                                                       glucocorticoids+
                                                                       plan_class+
                                                                       index_period'),
                                           dataset=rslt$cov_vs_oriac,
                                           description='iptw_covid_oriac') {
  
  dataset_final <- set_iptw_levels(iptw_data=dataset)
  
  #' Propensity Score Model
  ps_model <- gbm(weight_formula,
                  distribution="bernoulli",
                  n.trees = 100,
                  interaction.depth = 4,
                  train.fraction = .8,
                  shrinkage = .0005,
                  data=dataset_final)
  
  pscore <- ps_model$fitted.values
  
  dataset_final$p_score <- predict(ps_model, type = "response")
  
  visuals = list()
  viz_pscores <- dataset_final %>% 
    select(person_id, sub_cohort, covid, p_score) %>% 
    ggplot() + 
    geom_boxplot(aes(x=sub_cohort, y=p_score)) +
    theme_bw() +
    scale_color_brewer(palette = "Dark2") +
    theme(legend.position = "top")
  
  visuals$uncapped_pscores <- viz_pscores
  
  lower_cut <- dataset_final %>% filter(covid==1L) %>% pull(p_score) %>% quantile(0.005)
  upper_cut <- dataset_final %>% filter(covid==0L) %>% pull(p_score) %>% quantile(0.995)
  
  dataset_final_trim <- dataset_final %>%
    filter((p_score > lower_cut) & (p_score < upper_cut))
  
  ps_model_trim <- glm(weight_formula,
                       family=binomial,
                       data=dataset_final_trim)
  
  pscore2 <- ps_model_trim$fitted.values
  
  dataset_final_trim$p_score <- predict(ps_model_trim, type = "response")
  
  viz_pscores_trim <- dataset_final_trim %>% 
    select(person_id, sub_cohort, covid, p_score) %>% 
    ggplot() + 
    geom_boxplot(aes(x=sub_cohort, y=p_score)) +
    theme_bw() +
    scale_color_brewer(palette = "Dark2") +
    theme(legend.position = "top")
  
  visuals$capped_pscores <- viz_pscores_trim
  
  
  #' Visualize overlap of raw propensity scores
  viz_ps_overlap <- ggplot(data=dataset_final_trim, aes(x=p_score, group=sub_cohort, fill=sub_cohort)) +
    geom_density(adjust=1.5, alpha=.4) +
    theme_bw() +
    scale_fill_brewer(palette = "Dark2") +
    theme(legend.position = "top")
  
  visuals$trimmed_ps_overlap <- viz_ps_overlap
  
  
  #' Insert weights to get traditional IPTW for each participant (ATE)
  dataset_final_trim$capped_iptw = ifelse((dataset_final_trim$covid==1L),
                                          (1/dataset_final_trim$p_score),
                                          (1/(1-dataset_final_trim$p_score))
  )
  
  ps_model_balanced <- glm(weight_formula, family=quasibinomial(link="logit"),
                           data = dataset_final_trim,
                           weights = capped_iptw)
  # Now see the covariates are not statistically significant anymore
  
  viz_pscore_v_iptw <-
    dataset_final_trim %>% 
    select(person_id, sub_cohort, covid, p_score, capped_iptw) %>% 
    ggplot() + 
    geom_point(aes(x=p_score, y = capped_iptw, color=sub_cohort)) +
    theme_bw() +
    scale_color_brewer(palette = "Dark2") +
    theme(legend.position = "top")
  
  visuals$trimmed_pscore_v_iptw <- viz_pscore_v_iptw
  
  
  #' Make a vector of the most extreme weights
  extreme_weights <- c('covid propensity score pre-trimming',
                       'control propensity score pre-trimming',
                       'covid propensity score post-trimming',
                       'control propensity score post-trimming',
                       'covid IPTW weight post-trimming',
                       'control IPTW weight post-trimming')
  
  extreme1 <- dataset_final %>% filter(covid==1L) %>% summarise(max_p=max(p_score)) %>% pull()
  extreme2 <- dataset_final %>% filter(covid==1L) %>% summarise(min_p=min(p_score)) %>% pull()
  extreme3 <- dataset_final %>% filter(covid==0L) %>% summarise(max_p=max(p_score)) %>% pull()
  extreme4 <- dataset_final %>% filter(covid==0L) %>% summarise(min_p=min(p_score)) %>% pull()
  extreme5 <- dataset_final_trim %>% filter(covid==1L) %>% summarise(max_p=max(p_score)) %>% pull()
  extreme6 <- dataset_final_trim %>% filter(covid==1L) %>% summarise(min_p=min(p_score)) %>% pull()
  extreme7 <- dataset_final_trim %>% filter(covid==0L) %>% summarise(max_p=max(p_score)) %>% pull()
  extreme8 <- dataset_final_trim %>% filter(covid==0L) %>% summarise(min_p=min(p_score)) %>% pull()
  extreme9 <- dataset_final_trim %>% filter(covid==1L) %>% summarise(max_iptw=max(capped_iptw)) %>% pull()
  extreme10 <- dataset_final_trim %>% filter(covid==1L) %>% summarise(min_iptw=min(capped_iptw)) %>% pull()
  extreme11 <- dataset_final_trim %>% filter(covid==0L) %>% summarise(max_iptw=max(capped_iptw)) %>% pull()
  extreme12 <- dataset_final_trim %>% filter(covid==0L) %>% summarise(min_iptw=min(capped_iptw)) %>% pull()
  extremes_high <- c(extreme1, extreme3, extreme5, extreme7, extreme9, extreme11)
  extremes_low <- c(extreme2, extreme4, extreme6, extreme8, extreme10, extreme12)
  
  extreme_final <- cbind(extreme_weights, extremes_high, extremes_low) %>%
    as_tibble() %>%
    mutate(extremes_high=as.numeric(extremes_high),
           extremes_low=as.numeric(extremes_low))
  
  #' Note that this function is producing the exact same weights as the manual insertion in the ifelse statement above
  #' (final_dataset_trimmed$capped_iptw and w=w.out1$weights are identical)
  #' But weightit has nice love plots and balance tables that can be used for visualization
  w.out1 <- WeightIt::weightit(weight_formula,
                               data = dataset_final_trim,
                               estimand = "ATE", method = "ps")
  
  set.cobalt.options(binary = "std")
  
  weights_auto <- summary(w.out1)
  
  final <- list()
  #final$ir_thou <- ir_thou
  final$dataset <- dataset_final_trim # gets the final dataset with the propensity scores using by-hand computation of the logistic regression
  
  final$love_plot <- love.plot(w.out1) +
    scale_color_brewer(palette = "Dark2") +
    theme(legend.position = "top")# gets the love plot weights using weightIt
  
  final$balance_tbl <- bal.tab(w.out1,
                               binary='std', continuous='std')
  final$smds <- final$balance_tbl$Balance %>% select(-Diff.Un) %>%
    rename(Diff_Adj=`Diff.Adj`) %>%
    mutate(SMD_gt_than_pt2=case_when(abs(Diff_Adj) > 0.2 ~ 1L,
                                     TRUE ~ 0L))
  final$extreme_final <- extreme_final
  
  
  
  # TODO find some way to save all the plots in one document or something
  visuals_plot_sheet <- ggarrange(visuals$uncapped_pscores, visuals$capped_pscores,
                                  visuals$trimmed_ps_overlap, visuals$trimmed_pscore_v_iptw,
                                  final$love_plot,
                                  ncol = 2, nrow = 3, labels = "AUTO")
  ## Save the plot outcome as a summary sheet in results 
  ggsave(paste0("results/visual_summary_sheet_iptw_", description, ".png"), 
         plot = visuals_plot_sheet,
         width = 8.5, height = 11, dpi = 500
  )
  ggsave(paste0("results/visual_summary_sheet_iptw_", description, ".svg"), 
         plot = visuals_plot_sheet,
         width = 8.5, height = 11, dpi = 500
  )
  ggsave(paste0("results/visual_summary_sheet_iptw_", description, ".pdf"), 
         plot = visuals_plot_sheet,
         width = 8.5, height = 11, dpi = 500
  )
  
  #' Output final IPTW table with weights, post trimming
  output_tbl(dataset_final_trim, paste0(description,'_tabc'), indexes=c('person_id'))
  
  return(final)
  
}






#' Use WeightIt package to make weights for IPTW; and add MW
#' @param weight_formula formula to weight your model; make sure outcome is entered as binary 1/0
#' @param dataset dataset for weighting
#' @param description string descriptor for output
#' @return dataset with propensity scores, IPTW weights, IPTW & matching weights (MW), and 95th percentile-capped IPTW weights
#' samplesize (sample size), weightrange (highest and lowest weights), weighttop (highest weights),
#' coefofvar (coefficient of variation), iptw_blance (balance table for iptw), iptw_observations (iptw obserations),
#' iptw (dataset with capped weights?) maxweights_postcap (maximum weights post 95% weight capping)
make_iptw_and_mw_weights <- function(weight_formula=as.formula('covid~age_group+
                                                                       sex_cat+
                                                                       raceth_cat+
                                                                       total_pre+
                                                                       complex_chronic_flag+
                                                                       prior_sud+
                                                                       bmi_cat+
                                                                       glucocorticoids+
                                                                       plan_class+
                                                                       index_period'),
                                     dataset=rslt$cov_vs_oriac,
                                     description='iptw_covid_oriac') {
  
  library(WeightIt)
  library(cobalt)
  
  dataset_final <- dataset #%>% filter(site %in% site_list)
  
  #' pre-weighting SMDs
  
  balance_table_pre <- bal.tab(weight_formula, data=dataset_final,
                               estimand = 'ATT',
                               m.threshold=.1, disp.v.ratio=TRUE)
  
  balance_info_pre <-
    balance_table_pre$Balance %>% 
    rownames_to_column('variable') %>%
    as_tibble() %>%
    output_tbl(paste0(description,'_balance_pre'))
  
  #' weighting
  w.out <- weightit(weight_formula,
                    data=dataset_final, verbose=TRUE, method="glm", estimand="ATE")
  
  full_cohort_iptw <- dataset_final
  full_cohort_iptw$iptw_weights <- w.out$weights
  full_cohort_iptw$ps <- w.out$ps
  
  #' Modify IPTW to MW form (multiply numerator by min(1-ps_i, ps_i))
  full_cohort_iptw$iptw_mw <- full_cohort_iptw$iptw_weights * pmin(1-full_cohort_iptw$ps, full_cohort_iptw$ps)
  
  #' plot raw propensity scores
  # ggplot(data=full_cohort_iptw, aes(x=ps, group=covid, fill=covid)) +
  #   geom_density(adjust=1.5, alpha=0.4)
  
  #' summary:
  w.out_sum<-summary(w.out)
  
  #sample_size <- w.out_sum$effective.sample_size %>% as_tibble() %>% output_tbl(paste0(description,'_samplesize'))
  weight_range <- w.out_sum$weight.range %>% as_tibble() %>% output_tbl(paste0(description,'_weightrange'))
  biggest_weights_pre_capping <- w.out_sum$weight.top %>% as_tibble() %>% output_tbl(paste0(description,'_weighttop'))
  coef_of_var <- w.out_sum$coef.of.var %>% as_tibble() %>% output_tbl(paste0(description,'_coefofvar'))
  
  #' balance table:
  balance_table <- bal.tab(w.out, m.threshold=.1, disp.v.ratio=TRUE)
  
  balance_info <-
    balance_table$Balance %>% 
    rownames_to_column('variable') %>%
    as_tibble() %>%
    output_tbl(paste0(description,'_balance'))
  
  balance_table$Observations %>% output_tbl(paste0(description,'_observations'))
  
  #' 95th percentile capping on iptw weights only:
  #' cap weights to 95th percentile
  cutoff_upper <- quantile(full_cohort_iptw$iptw_weights, 0.995)
  cutoff_upper_num <- cutoff_upper[[1]]
  
  cutoff_lower <- quantile(full_cohort_iptw$iptw_weights, 0.005)
  cutoff_lower_num <- cutoff_lower[[1]]
  
  #' adding weights
  full_cohort_iptw_final <- full_cohort_iptw %>%
    mutate(capped_iptw=case_when(
      iptw_weights>=cutoff_upper_num~cutoff_upper_num,
      iptw_weights<=cutoff_lower_num~cutoff_lower_num,
      TRUE~iptw_weights
    ))
  
  output_tbl(full_cohort_iptw_final,
             paste0(description,'_table'),
             indexes=c('person_id'))
  
  min_weights <- full_cohort_iptw_final %>%
    group_by(covid) %>%
    summarise(min_iptw_weight=min(iptw_weights),
              min_capped_weight=min(capped_iptw)) %>%
    output_tbl(paste0(description,'_minweights_postcap'))
  
  max_weights <- full_cohort_iptw_final %>%
    group_by(covid) %>%
    summarise(max_iptw_weight=max(iptw_weights),
              max_capped_weight=max(capped_iptw)) %>%
    output_tbl(paste0(description,'_maxweights_postcap'))
  
  return(full_cohort_iptw_final)
  
}


#' Set reference levels for IPTW models
set_iptw_levels <- function(iptw_data=rslt$iptw_data) {
  
  iptw_data$age_group <- factor(iptw_data$age_cat_at_ce) %>% relevel(ref = '4-5y old')
  iptw_data$sex_cat <- factor(iptw_data$sex_cat) %>% relevel(ref = 'Male')
  iptw_data$race_eth_cat <- factor(iptw_data$race_eth_cat) %>% relevel(ref = 'NH_White')
  iptw_data$util_other <- factor(iptw_data$util_other) %>% relevel(ref='low_utilizer')
  iptw_data$util_inpatient <- factor(iptw_data$util_inpatient) %>% relevel(ref='no_visits')
  iptw_data$util_outpatient <- factor(iptw_data$util_outpatient) %>% relevel(ref='low_utilizer')
  iptw_data$util_ed <- factor(iptw_data$util_ed) %>% relevel(ref='no_visits')
  #iptw_data$total_pre <- factor(iptw_data$total_pre) %>% relevel(ref = '01 to 05 visits')
  # iptw_data$total_pre_1y <- factor(iptw_data$total_pre_1y) %>% relevel(ref = '01_to_04_visits')
  #iptw_data$visit_type <- factor(iptw_data$visit_type) %>% relevel(ref = 'Outpatient_and_telemedicine')
  # iptw_data$index_period2 <- factor(iptw_data$index_period2) %>% relevel(ref = '01_Mar_Aug_2020')
  # iptw_data$complex_chronic_flag <- factor(iptw_data$complex_chronic_flag) %>% relevel(ref = 'not_chronic_or_complex_chronic')
  #iptw_data$plan_class <- factor(iptw_data$plan_class) %>% relevel(ref = 'Private')
  
  return(iptw_data)
}





#' Get iptw tables for secondary outcomes (abx_iv, ed_visits, imaging, suppox, hba1c)
#' Weights by a count of occurrences in the prior 6 months
get_all_iptw_secondary_tables <- function(scd_subtypes=c('HbSS_and_HbS_beta_zero','HbSC','Other_Nonspecific'),
                                          subtype_str='all',
                                          exclude_sites=c('mshs'),
                                          outcome_str_list=c('abx_iv','ed_visits','imaging','suppox','hba1c'),
                                          cohort_tbl=results_tbl('scdlong_wcovs_oriac_outcomes')) {
  
  #' Loop through list of outcome strings
  for (i in 1:length(outcome_str_list)) {
    curr_outcome <- outcome_str_list[[i]]
    
    #' Set loop k to run through covid-negative control comparisons
    ctrl_group <- list('scd_covid_neg_test',
                       'scd_covid_neg_noev',
                       'scd_covid_neg_resp')
    
    #' Set loop j to run through list of different time periods (1y, pd1, pd2)
    period_list_min <- list(28L, 28L, 180L)
    period_list_max <- list(365L, 179L, 365L)
    date_str <- list('28_to_365','28_to_179','180_to_365')
    
    period_list_desc <- list('1y',
                             'pd1',
                             'pd2')
    
    #' Make name of outcome in 6 months prior to index date
    prev_info <- paste0(outcome_str_list[[i]],'_pre_179_to_28_ct')
    
    #' Loop through covid-negative control comparisons
    for (k in 1:length(cohort_list)) {
      
      ctrl_group_desc <- gsub("scd_covid_","",ctrl_group[[k]]) #' description of control group
      
      curr_ctrl_group <- ctrl_group[[k]]
      
      #' Ready data to run for iptw weights
      cohort_final <- cohort_tbl %>%
        filter(!site %in% exclude_sites) %>%
        filter(group %in% c('scd_with_covid',curr_ctrl_group)) %>%
        mutate(covid=case_when(group=='scd_with_covid' ~ 1L,
                               group==curr_ctrl_group ~ 0L)) %>%
        filter(scd_type %in% scd_subtypes) %>%
        collect()
      
      weight_formula <- as.formula(paste0('covid~age_group+
                               sex_cat+
                               raceth_cat+
                               total_pre+
                               complex_chronic_flag+
                               index_period2+
                               site+',prev_info))
      
      desc <- paste0('iptw_',ctrl_group_desc,'_',subtype_str,'_',curr_outcome)
      
      iptw_covid_neg_noev_other <- make_iptw_and_mw_weights_colby(weight_formula=weight_formula,
                                                                  dataset=cohort_final,
                                                                  description=desc)
      
      print(desc)
      
    }
    
  }
  
}




run_logistic_regression <- function(df,
                                    model_formula,
                                    weight_method = "none") {
  if (weight_method=="lr") {
    print("using LR weights")
    model_fit <- glm(formula =as.formula(model_formula),
                     data = df, 
                     family = quasibinomial(link="logit"),
                     weights = capped_iptw) 
  } else if (weight_method=="gbm") {
    model_fit <- glm(formula = as.formula(model_formula),
                     data = df, 
                     family = quasibinomial(link="logit"),
                     weights = iptw_gbm)
  } else{
    model_fit <- glm(formula = as.formula(model_formula),
                     data = df, 
                     family = binomial(link="logit"))
  }
  
  model_fit %>% return()
}
