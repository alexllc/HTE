library("colorout")

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
source("./HTE_expression.R")

usrwd = "/exeh_4/alex_lau"
setwd(paste0(usrwd, "/HTE/wd/expression_HTE"))


perm_genes = c(
# "ENSG00000245149","ENSG00000151012","ENSG00000153233","ENSG00000251692","ENSG00000169218","ENSG00000182557",
"ENSG00000138379","ENSG00000186204","ENSG00000105278","ENSG00000117148","ENSG00000188176","ENSG00000112214","ENSG00000116544","ENSG00000120049","ENSG00000267365","ENSG00000129749")

bins = seq(0,0.75,0.25)

for(gene in perm_genes) {

    print(paste0("treatment: ", gene))

    gene_ex = whole_dataset[,gene]

    for (i in bins){
        print(paste0("quantile: ", i))
        binary_tx = as.numeric(gene_ex > quantile(gene_ex, i) & gene_ex < quantile(gene_ex, i+0.25))
    # print(binary_tx)
        if (length(unique(binary_tx)) == 1 | sum(binary_tx) < length(binary_tx)*0.1) {
                        print("Gene expression distribution too sparse, skipping.")
                        next
                    }
        mod_ds = cbind(whole_dataset, binary_tx)
        mod_ds = dplyr::select(mod_ds, -gene)
        names(mod_ds)[names(mod_ds)=="binary_tx"] = gene
        mod_covar = dplyr::select(mod_ds, -c("donorId", "outcome"))
        obsNumber <- dim(covar_mat)[1]
        trainId <- sample(1: obsNumber, floor(obsNumber/2), replace = FALSE)
        registerDoParallel(10)

        output_file = paste0("./result/HNSC/alt_dirt/bin_test/", i, "_")

        result <- run.hte(mod_covar, gene, mod_ds, project, covar_type = "expression", trainId, seed = 111, is.binary = F, is_save = T, save_split = T, is.tuned = F, thres = 0.75, n_core = 8, output_directory = output_file)
        write.csv(result[[1]], paste0(output_file, project, '_expression_correlation_test_result.csv'), quote = F, row.names = F)
        write.csv(result[[2]], paste0(output_file, project, '_expression_calibration_result.csv'), quote = F, row.names = F)
        write.csv(result[[3]], paste0(output_file, project, '_expression_median_t_test_result.csv'), quote = F, row.names = F)
        write.csv(result[[4]], paste0(output_file, project, '_expression_permutate_testing_result.csv'), quote = F, row.names = F)

    }

}