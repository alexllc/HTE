# 1. Please make sure you have installed all the packages below before running
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

# 2. Make sure all four accompanying scripts are in the same directory as the header script
source("./grf_parameters.R")
source("./HTE_main_functions.R")
source("./HTE_validation_functions.R")
source("./survival_imputation.R")

# 3. Set your project and make sure you have created a file called result under the working directory
project = NULL
output_file = paste0("./result/", project, "/")

# 4. Prepare covariate matrix, whole dataset matrix and a vectoor of treatment types
covar_mat = NULL
whole_dataset = NULL
tx_vector = NULL

obsNumber <- dim(covar_mat)[1]
trainId <- sample(1: obsNumber, floor(obsNumber/2), replace = FALSE)
registerDoParallel(30)

result <- run.hte(covariatcovar_mat, tx_veector, all_data, trainId, seed = 111, is.binary = T, is_save = T, save_split = T, is.tuned = F, thres = 0.75, n_core = 8, output_directory = output_file)
write.csv(result[[1]], paste0(output_file, '_expression_correlation_test_result.csv'), quote = F, row.names = F)
write.csv(result[[2]], paste0(output_file, '_expression_calibration_result.csv'), quote = F, row.names = F)
write.csv(result[[3]], paste0(output_file, '_expression_median_t_test_result.csv'), quote = F, row.names = F)
write.csv(result[[4]], paste0(output_file, '_expression_permutate_testing_result.csv'), quote = F, row.names = F)