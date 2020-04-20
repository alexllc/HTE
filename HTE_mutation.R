# 1. Please make sure you have installed all the packages below before running
# For main function
library(doParallel)
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
library(readxl)
library(survival)

# For expression retreival
library(TCGAbiolinks)
library(DT)
library(SummarizedExperiment)
library(plyr)
library(dplyr)
library(tidyr)
library(biomaRt)
library(RTCGAToolbox)

# 2. Make sure all four accompanying scripts are in the same directory as the header script
#setwd("./HTE")
source("./grf_parameters.R")
source("./HTE_main_functions.R")
source("./HTE_validation_functions.R")
source("./survival_imputation.R")

#usrwd = "/exeh_4/alex_lau"
usrwd = "/home/alex/project"
setwd(paste0(usrwd, "/HTE/wd/mut_HTE"))

project = "BRCA"

output_file = paste0("./result/", project, "/")


## Clinical data
if (!file.exists("../TCGA_CDR_clean.csv")) {

    download.file("https://ars.els-cdn.com/content/image/1-s2.0-S0092867418302290-mmc1.xlsx", "TCGA-CDR-SupplementalTableS1.xlsx")

    cdr = read_excel("TCGA-CDR-SupplementalTableS1.xlsx")
    clinical_dat = dplyr::select(cdr, c(bcr_patient_barcode, type, age_at_initial_pathologic_diagnosis,  gender, ajcc_pathologic_tumor_stage, OS, OS.time))

    #patient <- patient[patient$OS.time!=0,]
    clinical_dat[clinical_dat == "#N/A"] <- NA
    clinical_dat <- subset(clinical_dat, !is.na(OS.time) & !is.na(OS) & !is.na(age_at_initial_pathologic_diagnosis))
    colnames(clinical_dat)[colnames(clinical_dat)=="bcr_patient_barcode"] = "donorId"
# This renaming step is critical as the HTE main function will rely on the column named outcome to indicate Y
    #patient$gender = as.numeric(patient$gender == 'FEMALE')

    write.csv(clinical_dat, "TCGA_CDR_clean.csv", row.names = F)
} else {
    clinical_dat = read.csv("TCGA_CDR_clean.csv")
}

# specify cancer type here
ss_patient <- subset(clinical_dat, type %in% project & OS.time > 0)
surv.times <- as.numeric(ss_patient$OS.time)
cens <- as.numeric(ss_patient$OS)

# get imputed log survival times
max.censored <- max(surv.times[cens == 0])
cens[surv.times == max.censored] <- 1
outcome = impute.survival(surv.times, cens)

# attach imputed.log.times to original dataset
ss_patient <- cbind(ss_patient, outcome)

for (c in colnames(ss_patient)) {

    if (!is.numeric(ss_patient[,c]) && c != "donorId") {
        which.one <- which( levels(ss_patient[,c]) == "")
        levels(ss_patient[,c])[which.one] <- NA
        ss_patient[,c] = sapply(sapply(ss_patient[,c], as.factor), as.numeric) 
        print(paste0(c, " is altered")) 
    }
}
ss_patient = dplyr::select(ss_patient, -c(type, OS, OS.time))
print("Processed patient dataframe: ")
head(ss_patient)

## mutation data

Pns_mat = read.csv(paste0("./Pns/TCGA-", project, "_Pns.csv"))

whole_dataset = left_join(ss_patient, Pns_mat, by = "donorId")
whole_dataset = whole_dataset[complete.cases(whole_dataset),]
covar_mat= dplyr::select(whole_dataset, -c("donorId", "outcome"))

tx_vector = colnames(covar_mat[,4:ncol(covar_mat)])
check = unlist(lapply(covar_mat[,4:ncol(covar_mat)], function(x) (length(unique(x)) == 1 | sum(x != 0) < length(x)*0.01))) # onyl use genes with at least 1% pop has the mutaion
tx_vector = tx_vector[!as.logical(check)]

obsNumber <- dim(covar_mat)[1]
trainId <- sample(1: obsNumber, floor(obsNumber/2), replace = FALSE)
registerDoParallel(10)

result <- run.hte(covar_mat, tx_vector, whole_dataset, project, covar_type = "mutation", trainId, seed = 111, is.binary = T, is_save = T, save_split = T, is.tuned = F, thres = 0.75, n_core = 8, output_directory = output_file)
write.csv(result[[1]], paste0(output_file, project, '_expression_correlation_test_result.csv'), quote = F, row.names = F)
write.csv(result[[2]], paste0(output_file, project, '_expression_calibration_result.csv'), quote = F, row.names = F)
write.csv(result[[3]], paste0(output_file, project, '_expression_median_t_test_result.csv'), quote = F, row.names = F)
write.csv(result[[4]], paste0(output_file, project, '_expression_permutate_testing_result.csv'), quote = F, row.names = F)