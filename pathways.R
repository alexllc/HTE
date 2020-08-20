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


format_tcga_patient <- function(pat_ls) {
    tmp = strsplit(pat_ls, "-")
    tmp = unlist(lapply(tmp, function(x) paste(x[[1]], x[[2]], x[[3]], sep = "-")))
    return(unlist(tmp))
}

props = read.table("/exeh_4/alex_lau/proj/HTE/wd/pathways/props_pathway_score.tsv", header = T, sep = '\t')
tmp = strsplit(colnames(props), '\\.')
tmp = sapply(tmp, function(x) x[[1]])
colnames(props) = tmp
props$donorId = format_tcga_patient(rownames(props))


cdr = read.csv("/exeh_4/alex_lau/proj/HTE/wd/TCGA_CDR_clean.csv")
cdr = dplyr::filter(cdr, type == "BRCA")
labels = c("[Discrepancy]","[Not Applicable]","[Not Available]","[Unknown]")
cdr$ajcc_pathologic_tumor_stage[which(cdr$ajcc_pathologic_tumor_stage %in% labels)] = NA

for (c in colnames(cdr)) {
    if (!is.numeric(cdr[,c]) && c != "donorId") {
        which.one <- which( levels(cdr[,c]) == "")
        levels(cdr[,c])[which.one] <- NA
        cdr[,c] = sapply(sapply(cdr[,c], as.factor), as.numeric) 
        print(paste0(c, " is altered")) 
    }
}
attach(cdr)
tcga_imp = NNMIS(ajcc_pathologic_tumor_stage, xa = age_at_initial_pathologic_diagnosis, xb = age_at_initial_pathologic_diagnosis, time = OS.time, event = OS, imputeCT = T, Seed = 2020, mc.cores = 60)
detach(cdr)
tcga_imp_surv = tcga_imp$dat.T.NNMI %>% mutate(mean = rowMeans(.))
cdr$outcome = tcga_imp_surv$mean

whole_dat = left_join(cdr, props, by = "donorId") %>% dplyr::select(-c(type, OS, OS.time))
whole_dat = whole_dat[complete.cases(whole_dat),]
covar_mat = dplyr::select(whole_dat, -c(donorId, outcome))


obsNumber <- dim(covar_mat)[1]
trainId <- sample(1: obsNumber, floor(obsNumber/2), replace = FALSE)
registerDoParallel(10)
output_file = "/exeh_4/alex_lau/proj/HTE/wd/pathways/results/PROPS_"
project = "BRCA"

result <- run.hte(covar_mat, tx_vector, whole_dat, project, covar_type = "UQ", txdirct = NULL, trainId, seed = 111, is.binary = F, is_save = T, save_split = T, is.tuned = F, thres = 0.75, n_core = 8, output_directory = output_file)
write.csv(result[[1]], paste0(output_file, project, '_expression_correlation_test_result.csv'), quote = F, row.names = F)
write.csv(result[[2]], paste0(output_file, project, '_expression_calibration_result.csv'), quote = F, row.names = F)
write.csv(result[[3]], paste0(output_file, project, '_expression_median_t_test_result.csv'), quote = F, row.names = F)
write.csv(result[[4]], paste0(output_file, project, '_expression_permutate_testing_result.csv'), quote = F, row.names = F)
