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
# library(TCGAbiolinks)
library(DT)
library(SummarizedExperiment)
library(plyr)
library(dplyr)
library(tidyr)
library(biomaRt)
# library(RTCGAToolbox)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(Homo.sapiens)

# 2. Make sure all four accompanying scripts are in the same directory as the header script
# setwd("../../HTE")
source("./grf_parameters.R")
source("./HTE_main_functions.R")
source("./HTE_validation_functions.R")
source("./survival_imputation.R")
# setwd("../expression_HTE/METABRIC")
setwd("../mut_HTE/METABRIC")
# Prepare clinical data

# dirct = "./data_jul4/"
dirct = "./data_jun11/"
project = "BRCA"
output_file = paste0("./result/", project, "/")

samp = fread(paste0(dirct, 'data_clinical_sample.txt'), skip=4)
#  [1] "#Identifier to uniquely specify a patient."
#  [2] "A unique sample identifier." 
#  [3] "Cancer Type" 
#  [4] "Cancer Type Detailed"
#  [5] "ER Status" 
#  [6] "HER2 Status" 
#  [7] "Numeric value to express the degree of abnormality of cancer cells, a measure of differentiation and aggressiveness."
#  [8] "Oncotree Code" 
#  [9] "PR Status" 
# [10] "The type of sample (i.e., normal, primary, met, recurrence)."
# [11] "Tumor size." 
# [12] "Tumor stage."
samp = samp[,c(1,4:7,9)]

pat = fread(paste0(dirct, 'data_clinical_patient.txt'), skip=4)
#[1] "#Identifier to uniquely specify a patient. 
#[2] "Number of lymphnodes positive"
#[3] "Nottingham prognostic index"
#[4] "Tumor Content"
#[5] "Chemotherapy."
#[6] "Cohort."
#[7] "ER status measured by IHC"
#[8] "HER2 status measured by SNP6" 
#[9] "Hormone Therapy"
# [10] "Inferred Menopausal State"
# [11] "Integrative Cluster"
# [12] "Age at Diagnosis" 
# [13] "Overall survival in months since initial diagonosis." 
# [14] "Overall patient survival status." 
# [15] "Pam50 + Claudin-low subtype"
# [16] "3-Gene classifier subtype"
# [17] "The survival state of the person."
# [18] "For tumors in paired organs, designates the side on which the cancer originates." 
# [19] "Radio Therapy"
# [20] "Text to describe a tumor's histologic subtype or mixed diagnosis that is different from previously specified options."
# [21] "Type of Breast Surgery" (almost all get masectomy)

pat = pat[,c(1,2,5,7:10,12,13,14,19)]
pat = pat[which(!(pat$OS_MONTHS <= 0 | is.na(pat$OS_MONTHS))),]
pat$OS_STATUS = as.integer(pat$OS_STATUS=="1:DECEASED") # METABRIC updated the OS status term
# Convert OS mo to days
pat$OS_MONTHS = pat$OS_MONTHS * 30.4167
# Impute survival time
max_censored = max(pat$OS_MONTHS[pat$OS_STATUS == 0], na.rm = T)
pat$OS_STATUS[pat$OS_MONTHS==max_censored] = 1
pat$outcome = impute.survival(pat$OS_MONTHS,  pat$OS_STATUS)
pat$OS_MONTHS = NULL
pat$OS_STATUS = NULL

# get imputed log survival times
clinical = merge(pat, samp, by = "PATIENT_ID")
clinical = as.data.frame(clinical[complete.cases(clinical),])

# Convert categorical variables into numeric
for (c in colnames(clinical)) {

    if (!is.numeric(clinical[,c]) && c != "PATIENT_ID") {
        which.one <- which( levels(clinical[,c]) == "")
        levels(clinical[,c])[which.one] <- NA
        clinical[,c] = sapply(sapply(clinical[,c], as.factor), as.numeric) 
        print(paste0(c, " is converted to numeric.")) 
    }
}

# Prepare microarray data
if(!file.exists("./transposed_data_expression_median.txt")) {

    medexp = fread(paste0(dirct,"data_expression_median.txt"))
    message("Transposing expression matrix, this could take a *LONG* while.")
    medexp = t(medexp)
    colnames(medexp) = medexp[1,]
    medexp = medexp[-c(1:2),]
    class(medexp) = "numeric"
    medexp = as.data.frame(medexp)
    medexp$PATIENT_ID = rownames(medexp)
    write.table(medexp, "transposed_data_expression_median.txt", sep = "\t", row.names = F) # as of R4.0 R no longer automatically converts strings as factors
} else {
    medexp = fread("transposed_data_expression_median.txt")
}


# HTE datset
wholedat = left_join(clinical, medexp, by = "PATIENT_ID")
colnames(wholedat)[1] = "donorId"
wholedat = wholedat[complete.cases(wholedat),]
covar = dplyr::select(wholedat, -c(donorId, outcome))


# add tx direction
cna = fread(paste0(dirct, "data_CNA.txt"))
cna = as.data.frame(t(as.matrix(cna)))
colnames(cna) = cna[1,]
cna = cna[-c(1,2),]
cna = cna[complete.cases(cna),]
cna = sapply(cna, as.numeric)
DEGs = apply(cna, 2, sum)


# retreive TCGA results
TCGA_cor_res = read.csv("/home/alex/project/HTE/wd/expression_HTE/result/DEAcor/BRCA/BRCA_expression_correlation_test_result.csv")
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
txnames = AnnotationDbi::select(Homo.sapiens, keys = unique(TCGA_cor_res$gene), columns = "SYMBOL", keytype = "ENSEMBL", multiVals = "CharacterList")
TCGA_genes = left_join(TCGA_cor_res, txnames, by = c("gene" = "ENSEMBL"))
TCGA_genes = unique(TCGA_genes$SYMBOL)

tx_vector = colnames(medexp[,-ncol(medexp)])
tx_vector = tx_vector[tx_vector %in% colnames(cna) & tx_vector %in% TCGA_genes]
# test code
tx_vector = tx_vector[1:5]

obsNumber <- dim(covar)[1]
trainId <- sample(1: obsNumber, floor(obsNumber/2), replace = FALSE)
registerDoParallel(10)


result <- run.hte(covar, tx_vector, wholedat, project, covar_type = "expression", txdirct = DEGs, trainId, seed = 111, is.binary = T, is_save = T, save_split = T, is.tuned = F, thres = 0.75, n_core = 6, output_directory = output_file)
write.csv(result[[1]], paste0(output_file, project, '_expression_correlation_test_result.csv'), quote = F, row.names = F)
write.csv(result[[2]], paste0(output_file, project, '_expression_calibration_result.csv'), quote = F, row.names = F)
write.csv(result[[3]], paste0(output_file, project, '_expression_median_t_test_result.csv'), quote = F, row.names = F)
write.csv(result[[4]], paste0(output_file, project, '_expression_permutate_testing_result.csv'), quote = F, row.names = F)