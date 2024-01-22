# cov-subsequent-infections
Study focusing on whether patients with SARS-CoV-2 are more or less at risk for future infections of any kind compared with sets of control patients.

For latest study progress, see [www.google.com](`reporting/study_progress_live.html`) for an update.

## `code/` files 

All files directly in the `code/` folder marked with the prefix `cohort_...` are files with relevant functions to the development of the base cohort or analysis. The following sub folders contain driver files responsible for executing to code to create cohorts, perform EDA, create the analytic dataset and perform analysis.

### Exploration

Based on proposed cohort inclusion criteria, files in here support exploration of the data

Also include preliminary data quality and availability analyses


### Cohort

Files in here support the creation of the cohort for the study

### Codeset

Files in here support the creation and output of codesets to the `specs/` folder


### Variables

Files in this folder support the development of logic for variable extraction for this specific study from data.
Outputs may be in the form of tables output to the results schema.


### Analysis

Files in this folder support the analysis of the completed dataset for the study.
Complex analysis outputs will include tables output to the results schema.

This folder will also include code supporting the creation of figures and tables, though code to generate such figures and tables can be found in the `reporting/` folder.

## `reporting/` files



