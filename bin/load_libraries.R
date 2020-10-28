## Script to load all libraries in one place

# For TCGA_data_functions
library(readxl)
library(TCGAbiolinks)
library(dplyr)
library(tidyr)
library(NNMIS)

# For survival_imputation, not necessary if you use NNMIS, NNMIS will load this package on its own
library(survival)

# For HTE_main_function
library(grf)
library(MASS)
library(doMC)
library(data.table)
library(survminer)
library(doParallel)
library(methods)

# For validation_functions
library(DescTools)

# For gfr_parameters
library(grf)
library(BART)
library(ranger)
library(randomForestSRC)
library(randomForest)

# For overlap_of_varimp, requires R ver 4.0 or above, if you do not need to use varimp overlap assessments, skip this.
library(foreach)
library(TFisher) 
library(metap)
library(parallel)


# For METABRIC_external_validation_functions.R
library()