library(dplyr)
library(NNMIS)

library(grf)
library(MASS)
library(doMC)
library(data.table)
library(survminer)
library(doParallel)
library(methods)
library(readxl)

# For causal forest
library(grf)
library(BART)
library(ranger)
library(randomForestSRC)
library(randomForest)

# For survival imputation
library(survival)
library(NNMIS)

source("./grf_parameters.R")
source("./HTE_main_functions.R")
source("./HTE_validation_functions.R")
source("./NNMIS_survival_imputation.R")


format_tcga_patient <- function(pat_ls) {
    tmp = strsplit(pat_ls, "-")
    tmp = unlist(lapply(tmp, function(x) paste(x[[1]], x[[2]], x[[3]], sep = "-")))
    return(unlist(tmp))
}

cancer_list = c(
                'BLCA',
                'COAD',
                'BRCA',
                # 'LGG', # no NT
                'GBM',
                'STAD',
                'HNSC',
                'KIRC',
                'LUAD',
                'LUSC',
                #'OV', # no NT
                'PRAD',
                #'SKCM', # no NT
                'THCA',
                'UCEC',
                'ESCA')

source("./NNMIS_survival_imputation.R")
setwd("/exeh_4/alex_lau/proj/HTE/wd/pathways")

for (c in cancer_list) {

    props = read.csv(paste0("/exeh_4/alex_lau/proj/HTE/wd/pathways/path_scores/props_pathway_score_", c, ".csv"))
    colnames(props)[colnames(props) == "X"] = "donorId"
    tx_vector = colnames(props)[2:length(props)]

    cdr = read_excel("/exeh_4/alex_lau/proj/HTE/wd/TCGA-CDR-SupplementalTableS1.xlsx")
    cdr = dplyr::filter(cdr, type == c)
    clinical = dplyr::select(cdr, all_of(c("bcr_patient_barcode", "age_at_initial_pathologic_diagnosis", "gender", "OS.time", "OS", "ajcc_pathologic_tumor_stage", "tumor_status"))) %>% filter(!is.na(OS.time)) %>% as.data.frame()
    colnames(clinical)[colnames(clinical) == "bcr_patient_barcode"] = "donorId"
    clinical = impute_with_NNMIS(clinical)

    whole_dat = left_join(clinical, props, by = "donorId") %>% dplyr::select(-c(OS, OS.time))
    whole_dat = whole_dat[complete.cases(whole_dat),]
    covar_mat = dplyr::select(whole_dat, -c(donorId, outcome))


    obsNumber <- dim(covar_mat)[1]
    trainId <- sample(1: obsNumber, floor(obsNumber/2), replace = FALSE)
    registerDoParallel(10)
    output_file = paste0("/exeh_4/alex_lau/proj/HTE/wd/pathways/pancan_results/LQ_tx/", c, "/PROPS_")

    result <- run.hte(covar_mat, tx_vector, whole_dat, c, covar_type = "LQ", txdirct = NULL, trainId, seed = 111, is.binary = F, is_save = T, save_split = T, is.tuned = F, thres = 0.75, n_core = 8, output_directory = output_file)
    write.csv(result[[1]], paste0(output_file, project, '_expression_correlation_test_result.csv'), quote = F, row.names = F)
    write.csv(result[[2]], paste0(output_file, project, '_expression_calibration_result.csv'), quote = F, row.names = F)
    write.csv(result[[3]], paste0(output_file, project, '_expression_median_t_test_result.csv'), quote = F, row.names = F)
    write.csv(result[[4]], paste0(output_file, project, '_expression_permutate_testing_result.csv'), quote = F, row.names = F)
}