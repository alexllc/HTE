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
library(NNMIS)

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
#setwd("../HTE")
source("./grf_parameters.R")
source("./HTE_main_functions.R")
source("./HTE_validation_functions.R")
source("./NNMIS_survival_imputation.R")
# source("./survival_imputation.R")

# usrwd = "/exeh_4/alex_lau"
usrwd = "/home/alex/project"
setwd(paste0(usrwd, "/HTE/wd/mut_HTE"))

project = 'BRCA'
# for (project in cancer_list) {}
output_file = paste0("./result/", "BRCA_w_ICGC_small_forest/")


## Clinical data
if (!file.exists("../TCGA_CDR_clean.csv")) {

    download.file("https://ars.els-cdn.com/content/image/1-s2.0-S0092867418302290-mmc1.xlsx", "TCGA-CDR-SupplementalTableS1.xlsx")

    cdr = read_excel("TCGA-CDR-SupplementalTableS1.xlsx")
    clinical_dat = dplyr::select(cdr, c(bcr_patient_barcode, type, age_at_initial_pathologic_diagnosis,  gender, ajcc_pathologic_tumor_stage, tumor_status, OS, OS.time))

    #patient <- patient[patient$OS.time!=0,]
    clinical_dat[clinical_dat == "#N/A"] <- NA
    clinical_dat <- subset(clinical_dat, !is.na(OS.time) & !is.na(OS) & !is.na(age_at_initial_pathologic_diagnosis))
    colnames(clinical_dat)[colnames(clinical_dat)=="bcr_patient_barcode"] = "donorId"
    clinical_dat = as.data.frame(clinical_dat)
# This renaming step is critical as the HTE main function will rely on the column named outcome to indicate Y
    #patient$gender = as.numeric(patient$gender == 'FEMALE')

    write.csv(clinical_dat, "TCGA_CDR_clean.csv", row.names = F)
} else {
    clinical_dat = read.csv("../TCGA_CDR_clean.csv")
}

# specify cancer type here
ss_patient <- subset(clinical_dat, type %in% project & OS.time > 0)
ss_patient = impute_with_NNMIS(ss_patient)
ss_patient = dplyr::select(ss_patient, -c(OS, OS.time))
print("Processed patient dataframe: ")
head(ss_patient)

#######################################################################
## mutation data
#######################################################################

# maf <- GDCquery_Maf(project, pipelines = "muse")
if (!file.exists(paste0(project, "_maf.csv"))) {
    m_query <- GDCquery(project = paste0("TCGA-", project),
                    data.category = "Simple Nucleotide Variation",
                    legacy = FALSE,
                    workflow.type = "MuSE Variant Aggregation and Masking",
                    data.type = "Masked Somatic Mutation"
                )

    GDCdownload(m_query)
    maf = GDCprepare(m_query)
    write.csv(maf, paste0("maf", project, "_maf.csv"), row.names = F)

} else {maf = read.csv(paste0("maf",project, "_maf.csv"))}


if (!file.exists(paste0(project, "_all_mut_freq.csv"))) {
pmaf = dplyr::filter(maf, BIOTYPE == "protein_coding")
# Sep 7@Alex: synonymous mutations are also counted 
pmaf = dplyr::select(pmaf, c(Hugo_Symbol, Tumor_Sample_Barcode))
pmaf = pmaf %>% group_by(Tumor_Sample_Barcode) %>% add_count(Hugo_Symbol)
pmaf = pmaf[!duplicated(pmaf),]
wtcga = pmaf %>% spread(Hugo_Symbol, n)
wtcga[is.na(wtcga)] = 0
wtcga$Tumor_Sample_Barcode = as.character(wtcga$Tumor_Sample_Barcode)
tmp = strsplit(wtcga$Tumor_Sample_Barcode, "-")
tmp = unlist(lapply(tmp, function(x) paste(x[[1]], x[[2]], x[[3]], sep = "-")))
wtcga$Tumor_Sample_Barcode <- unlist(tmp)
wtcga = as.data.frame(wtcga)
zero_counts = sapply(wtcga, function(x) sum(x == 0))
colnames(wtcga)[1] = "donorId"
write.csv(wtcga, paste0(project, "_all_mut_freq.csv"), row.names = F)
} else {
    wtcga = read.csv(paste0(project, "_all_mut_freq.csv"))
}


#### Incooperate ICGC data (or else it will be too sparse) Alex@sep14,2020
icgc = read.csv(paste0("./xena/ICGC_", project, "_df.csv"))
common_gene = intersect(colnames(icgc),colnames(wtcga))
icgc_sub = dplyr::select(icgc, all_of(c("sampleID", "donor_age_at_diagnosis", "donor_sex", "donor_tumour_stage_at_diagnosis", "outcome", "project_code", common_gene)))
tmp = strsplit(icgc_sub$project_code, "-")
tmp = sapply(tmp, function(x) x[[2]])
icgc_sub$project_code = tmp
icgc_sub$donor_sex = as.numeric(as.factor(icgc_sub$donor_sex))


tcga_sub = dplyr::select(wtcga, all_of(c("donorId", common_gene)))
ss_patient$country = "US"
tcga_sub = left_join(ss_patient, tcga_sub, by = "donorId")

common_clin = c("donorId", "age", "sex", "stage", "outcome", "country")
colnames(icgc_sub)[1:6] = common_clin
colnames(tcga_sub)[1:6] = common_clin

whole_dataset = rbind(icgc_sub, tcga_sub)
whole_dataset$country = as.numeric(as.factor(whole_dataset$country))
whole_dataset[is.na(whole_dataset)] = 0

# mskcc = read.table("./MSKCC-PRAD/data_mutations_extended.txt", sep  ='\t', header=T)
# smsk = dplyr::select(mskcc, Hugo_Symbol, Tumor_Sample_Barcode)
# smsk = smsk %>% group_by(Tumor_Sample_Barcode) %>% add_count(Hugo_Symbol)
# smsk = smsk[!duplicated(smsk),]
# wmsk = smsk %>% spread(Hugo_Symbol, n)
# wmsk[is.na(wmsk)] = 0

# cmon_gene = colnames(wmsk)[colnames(wmsk) %in% colnames(wtcga)]

# wtcga = dplyr::select(wtcga, c("donorId", cmon_gene))

# Pns_mat = read.csv(paste0("./Pns/TCGA-", project, "_Pns.csv"))


## Alex@sep14, 2020: no longer used
# whole_dataset = left_join(ss_patient, wtcga, by = "donorId")
# whole_dataset = whole_dataset[complete.cases(whole_dataset),]


freq_sum = sapply(whole_dataset[,7:dim(whole_dataset)[2]], function(x) sum(x != 0))
selection = freq_sum > dim(whole_dataset)[1]*0.02
whole_dataset = whole_dataset[,c(rep(TRUE, 6), selection)]
# whole_dataset = whole_dataset[complete.cases(whole_dataset),]
covar_mat= dplyr::select(whole_dataset, -c("donorId", "outcome"))
tx_vector = names(which(selection == T))

# check = unlist(lapply(covar_mat[,4:ncol(covar_mat)], function(x) (length(unique(x)) == 1 | sum(x != 0) < length(x)*0.01))) # onyl use genes with at least 1% pop has the mutaion
# tx_vector = tx_vector[!as.logical(check)]

obsNumber <- dim(covar_mat)[1]
trainId <- sample(1: obsNumber, floor(obsNumber/2), replace = FALSE)
registerDoParallel(10)

result <- run.hte(covar_mat, tx_vector, whole_dataset, project, covar_type = "mutation", trainId = trainId, seed = 111, is.binary = T, is_save = T, save_split = T, is.tuned = F, thres = 0.75, n_core = 8, output_directory = output_file, skip_perm = FALSE)
write.csv(result[[1]], paste0(output_file, project, '_expression_correlation_test_result.csv'), quote = F, row.names = F)
write.csv(result[[2]], paste0(output_file, project, '_expression_calibration_result.csv'), quote = F, row.names = F)
write.csv(result[[3]], paste0(output_file, project, '_expression_median_t_test_result.csv'), quote = F, row.names = F)
write.csv(result[[4]], paste0(output_file, project, '_expression_permutate_testing_result.csv'), quote = F, row.names = F)