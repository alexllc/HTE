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

source("./grf_parameters.R")
source("./HTE_main_functions.R")
source("./HTE_validation_functions.R")
source("./NNMIS_survival_imputation.R")
# source("./survival_imputation.R")

# usrwd = "/exeh_4/alex_lau"
usrwd = "/home/alex/project"
setwd(paste0(usrwd, "/HTE/wd/mut_HTE"))
project = "HNSC"

## Clinical data
if (!file.exists("../TCGA_CDR_clean.csv")) {

    download.file("https://ars.els-cdn.com/content/image/1-s2.0-S0092867418302290-mmc1.xlsx", "TCGA-CDR-SupplementalTableS1.xlsx")

    cdr = read_excel("TCGA-CDR-SupplementalTableS1.xlsx")
    clinical_dat = dplyr::select(cdr, c(bcr_patient_barcode, type, age_at_initial_pathologic_diagnosis,  gender, ajcc_pathologic_tumor_stage, tumor_status, OS, OS.time, tumor_status))

    clinical_dat[clinical_dat == "#N/A"] <- NA
    clinical_dat <- subset(clinical_dat, !is.na(OS.time) & !is.na(OS) & !is.na(age_at_initial_pathologic_diagnosis))
    colnames(clinical_dat)[colnames(clinical_dat)=="bcr_patient_barcode"] = "donorId"
    clinical_dat = as.data.frame(clinical_dat)
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
## mutation data import from TCGAbiolinks
#######################################################################

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

    # Selecting only protein coding genes & non-synonymous mutation for counts    
    pmaf = dplyr::filter(maf, BIOTYPE == "protein_coding")
    pmaf$SNNS = ifelse(pmaf$One_Consequence == "synonymous_variant" | is.na(pmaf$Amino_acids),"SN", "NS")
    NSpmaf = filter(pmaf, SNNS == "NS")
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

# Join two datasets
whole_dataset = left_join(ss_patient, wtcga, by = "donorId")
whole_dataset[is.na(whole_dataset)] = 0

freq_sum = sapply(whole_dataset[,7:dim(whole_dataset)[2]], function(x) sum(x != 0))
selection = freq_sum > dim(whole_dataset)[1]*0.02
whole_dataset = whole_dataset[,c(rep(TRUE, 6), selection)]
whole_dataset = whole_dataset[complete.cases(whole_dataset),]
covar_mat= dplyr::select(whole_dataset, -c("donorId", "outcome"))
tx_vector = names(which(selection == T))