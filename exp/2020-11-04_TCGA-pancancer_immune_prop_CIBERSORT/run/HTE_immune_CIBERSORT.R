#' Script to perform HTE using CIBERSORT cell types as treatment, expression as covaraites and OS as outcome
#' ***ATTENTION*** TO BE RUN AFTER grf_parameter.R and HTE_validation_functions.R scripts have been finalized
 
setwd("~/project/HTE/") # Run this script from the base directory

# source ./bin scripts
bin_ls = list.files("./bin")
for (bin in bin_ls){
    source(paste0("./bin/", bin))
}

# Specify cancer type
cancer_type = "HNSC"

# Fetch clinical data
clinical = fetch_clinical_data(cancer_type,)

# Fetch expression data
exp = fetch_exp_data(cancer_type, primaryTumorOnly = addBatch = TRUE, numericBatch = TRUE, scale = FALSE, primaryTumorOnly = FALSE, formatPatient = FALSE)
exp$donorId = rownames(exp)

# Load immune proportion dataframe as 'treatments'
############################# PLEASE COMPLETE ################################
proportion  = 

# Perform HTE
############################# PLEASE COMPLETE ################################

whole_dataset = 
covar_mat = 

obsNumber <- dim(covar_mat)[1]
trainId <- sample(1: obsNumber, floor(obsNumber/2), replace = FALSE)
registerDoParallel(10)

############################# PLEASE COMPLETE ################################
output_file = 

############################# PLEASE COMPLETE ################################
result <- run.hte(  )


write.csv(result[[1]], paste0(output_file, project, '_expression_correlation_test_result.csv'), quote = F, row.names = F)
write.csv(result[[2]], paste0(output_file, project, '_expression_calibration_result.csv'), quote = F, row.names = F)
write.csv(result[[3]], paste0(output_file, project, '_expression_median_t_test_result.csv'), quote = F, row.names = F)
write.csv(result[[4]], paste0(output_file, project, '_expression_permutate_testing_result.csv'), quote = F, row.names = F)