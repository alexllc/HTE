
# Setting up script for running CIBERSORT immune HTE

2020-11-04

***

## Script elements:

1. basic headers (set working directories and load libraries)
2. import clinical dataframe
3. import expression dataframe
4. import CIBERSORT immune proportion dataframe
5. call HTE function

## Corresponding update on files under ./bin

- `fetch_exp_data` function `@TCGA_data_functions.R` added option to keep all tissue types
- `fetch_clinical_data` function no longer keep only complete cases

## To be updated before running

(Larry) import CIBERSORT scores
(Larry) 