#' Script to run an HTE with mutation as treatment variable and expression as running covariates
# Run from base directory

# For data cleaning and processing
library(readxl)

# For main function
library(grf)
library(MASS)
library(doMC)
library(data.table)
library(survminer)
library(doParallel)
library(methods)

# For causal forest
library(grf)
library(BART)
library(ranger)
library(randomForestSRC)
library(randomForest)

# For survival imputation
library(survival)


# source bin scripts
bin_ls = list.files("./bin")
for (bin in bin_ls){
    source(paste0("./bin/", bin))
}

# Set param for this run
cancer_type = "BRCA"



# Retreive data
clinical = fetch_clinical_data(cancer_type)
mut = fetch_mut_data(cancer_type)
exp = fetch_exp_data(cancer_type, scale = FALSE)

# Construct HTE dataset

# Run HTE