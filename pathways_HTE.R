library(dplyr)
library(NNMIS)

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

for (c in cancer_list) {

    props = read.csv(paste0("/exeh_4/alex_lau/proj/HTE/wd/pathways/path_scores/props_pathway_score_", c, ".csv"))
    tmp = strsplit(colnames(props), '\\.')
    tmp = sapply(tmp, function(x) x[[1]])
    colnames(props) = tmp
    tx_vector = tmp
    props$donorId = format_tcga_patient(rownames(props))


    cdr = read.csv("/exeh_4/alex_lau/proj/HTE/wd/TCGA_CDR_clean.csv")
    cdr = dplyr::filter(cdr, type == c)
    

    whole_dat = left_join(cdr, props, by = "donorId") %>% dplyr::select(-c(type, OS, OS.time))
    whole_dat = whole_dat[complete.cases(whole_dat),]
    covar_mat = dplyr::select(whole_dat, -c(donorId, outcome))


    obsNumber <- dim(covar_mat)[1]
    trainId <- sample(1: obsNumber, floor(obsNumber/2), replace = FALSE)
    registerDoParallel(10)
    output_file = "/exeh_4/alex_lau/proj/HTE/wd/pathways/LQ_results/PROPS_"
    project = "BRCA"

    result <- run.hte(covar_mat, tx_vector, whole_dat, project, covar_type = "LQ", txdirct = NULL, trainId, seed = 111, is.binary = F, is_save = T, save_split = T, is.tuned = F, thres = 0.75, n_core = 8, output_directory = output_file)
    write.csv(result[[1]], paste0(output_file, project, '_expression_correlation_test_result.csv'), quote = F, row.names = F)
    write.csv(result[[2]], paste0(output_file, project, '_expression_calibration_result.csv'), quote = F, row.names = F)
    write.csv(result[[3]], paste0(output_file, project, '_expression_median_t_test_result.csv'), quote = F, row.names = F)
    write.csv(result[[4]], paste0(output_file, project, '_expression_permutate_testing_result.csv'), quote = F, row.names = F)
}