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
setwd(paste0(usrwd, "/HTE/wd/immune_HTE"))

project = HNSC
output_file = paste0("./result/", project, "_immune_out/")

## Prepare expression data
#######################################
### EXPRESSION DATA
#######################################
if (!file.exists(paste0("./HTSeqData/", project, "_exp.rda")) ) {

    # Fetch expression data from GDC
    g_query <- GDCquery(project = paste0("TCGA-", project),
                    data.category = "Transcriptome Profiling",
                    legacy = F,
                    data.type = "Gene Expression Quantification",
                    workflow.type = "HTSeq - FPKM-UQ",
                    sample.type = "Primary solid Tumor")

    GDCdownload(g_query)
    expdat <- GDCprepare(query = g_query,
                        save = TRUE,
                        save.filename = paste0("./HTSeqData/",project, "_exp.rda"))
    prep <- GDCprepare(g_query) 
} else {
    load(paste0("./HTSeqData/", project, "_exp.rda"))
    prep = data
}

exp_matrix <- SummarizedExperiment::assay(prep, "HTSeq - FPKM-UQ")

# Only select primary tumor samples
# Somehow this line removed so many patients
exp_matrix <- exp_matrix[,grep("-01.-", colnames(exp_matrix))]

#Scale FPKM value by patients
exp_matrix <- exp_matrix/apply(exp_matrix,2,max)
exp_matrix <- t(exp_matrix)

#Correct for batch effects
exp_matrix <- cbind(separate(as.data.frame(rownames(exp_matrix)), "rownames(exp_matrix)",
                        c(NA, "TSS", "patient", NA, "portion", "plate", "center"),  # skip var with NAs
                        sep = "-"), exp_matrix)
exp_matrix$bcr <- rownames(exp_matrix)

# Select only one aliquot
exp_matrix <- as.data.frame(exp_matrix %>% group_by(patient) %>% dplyr::slice(1))
rownames(exp_matrix) <- exp_matrix$bcr


batches = c("TSS", "patient", "portion", "plate", "center")
for (name in batches) {
    c = grep(name, colnames(exp_matrix))
    which.one <- which( levels(exp_matrix[,c]) == "")
    levels(exp_matrix[,c])[which.one] <- NA
    print(paste0(colnames(exp_matrix)[c], " is altered"))
    exp_matrix[,c] = sapply(sapply(exp_matrix[,c], as.factor), as.numeric) 
}

#Format patient
tmp = strsplit(rownames(exp_matrix), "-")
tmp = unlist(lapply(tmp, function(x) paste(x[[1]], x[[2]], x[[3]], sep = "-")))
rownames(exp_matrix) <- unlist(tmp)

exp_matrix$donorId <- rownames(exp_matrix)
exp_matrix <- dplyr::select(exp_matrix, -c(bcr, patient))


#######################################
### IMMUNE CELL PROPORTION DATA
#######################################
proportion = read.csv(paste0("./proportion/", project, "_immune_cells.csv"))

tx_vector = colnames(exp_matrix[5:ncol(exp_matrix)])

# 4. Prepare covariate matrix, whole dataset matrix and a vectoor of treatment types
whole_dataset = inner_join(exp_matrix, proportion, by = "donorId")
outcome = log(whole_dataset$CD8_Tcells)
#whole_dataset = dplyr::select(whole_dataset, c("donorId","outcome", "TSS", "portion", "plate", "center", tx_vector))
print("About to analysize the wholedataset:")
head(whole_dataset[,1:20])
covar_mat= dplyr::select(whole_dataset, -c("donorId", "CD8_Tcells"))

# 5. Run HTE with the whole dataset and covar matrix
obsNumber <- dim(covar_mat)[1]
trainId <- sample(1: obsNumber, floor(obsNumber/2), replace = FALSE)
registerDoParallel(10)

result <- run.hte(covar_mat, tx_vector, whole_dataset, project, covar_type = NULL, trainId, seed = 111, is.binary = F, is_save = T, save_split = T, is.tuned = F, thres = 0.75, n_core = 8, output_directory = output_file)
write.csv(result[[1]], paste0(output_file, project, '_expression_correlation_test_result.csv'), quote = F, row.names = F)
write.csv(result[[2]], paste0(output_file, project, '_expression_calibration_result.csv'), quote = F, row.names = F)
write.csv(result[[3]], paste0(output_file, project, '_expression_median_t_test_result.csv'), quote = F, row.names = F)
write.csv(result[[4]], paste0(output_file, project, '_expression_permutate_testing_result.csv'), quote = F, row.names = F)