# cov-subsequent-infections
Study focusing on whether patients with SARS-CoV-2 are more or less at risk for future infections of any kind compared with sets of control patients.

For latest study progress, see [`reporting/study_progress_live.html`](https://htmlpreview.github.io/?https://raw.githubusercontent.com/PEDSnet/cov-subsequent-infections/main/reporting/study_progress_live.html?token=GHSAT0AAAAAACNGFXCH4QURKPLRQUKLLO36ZNO6WFA) for an update.

## Driver files

Cohorts, covariates, and outcomes are generated in two main files:

1. `driver_rsv_study.R` contains the code for the study Aim 1: To assess whether patients with covid are more likely than patients with influenza or another respiratory infection to be at risk for subsequent RSV infection. In here, the cohort gets created, covariates computed, and outcomes generated, and finally an analytic dataset is created that can be seen from the end of the script.

2. `driver_respiratory_study.R` Similarly, this file contains the code for study Aim 2 and Aim 3. These aims share the same cohort, and so the file format is the same as the RSV study, with cohorts being generated, covariates computed, outcomes computed, and an analytic dataset created.

3. `driver_time_series_and_viz.R` contains some extra code to do two things: First, I did a computation of ALL the RSV cases during the 2022 RSV surge, just to get a sense of how big the surge was in total, and in order to support creating the figure for the landscape of infections. In this driver file you can find code to create that figure.

## Markdown and reporting files

Primary reporting files for current study analysis:

1. `rsv_progress_report.Rmd` contains the Table 1, attrition, and general information for the analysis of the Aim 1 RSV study. It includes descriptive statistics as well as sections on modeling and model comparison.

2. `resp_progress_report.Rmd` contains analogous tables and modeling for the Aim 2 and Aim 3 studies. 

## Table language

Some keywords to watch out for:

1. In general, keywords like "general infection" relate to Aim 3, "Respiratory" and "Resp" relate to Aim 2, and "RSV" to Aim 1.


### Exploration files

Archived code for study exploration.

Based on proposed cohort inclusion criteria, files in here support exploration of the data

Also include preliminary data quality and availability analyses


### Cohort

OLD/ARCHVIED: Files in here support the creation of the cohort for the study

### Codeset

Files in here support the creation and output of codesets to the `specs/` folder


### Variables

Files in this folder support the development of logic for variable extraction for this specific study from data.
Outputs may be in the form of tables output to the results schema.


### Analysis

OLD/ARCHVIED: Files in this folder support the analysis of the completed dataset for the study.
Complex analysis outputs will include tables output to the results schema.





