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

usrwd = "/exeh_4/alex_lau"
setwd(paste0(usrwd, "/HTE/wd/expression_HTE"))
# 3. Set your project and make sure you have created a file called result under the working directory
cancer_list = c(
                'BLCA',
                'COAD',
                'BRCA',
                #'LGG',
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

project = cancer_list[3]
output_file = paste0("./result/", project, "/")

# Fetch survival data from GDC
query <- GDCquery(project = paste0("TCGA-", project),
                data.category = "Clinical",
                file.type = "xml") # Obtain list from TCGAbiolinks:::getProjectSummary(<enter_project>)

#GDCdownload(query)

patient = GDCprepare_clinic(query, clinical.info = "patient")

save(patient, file = paste0("./clinical/", project, "_clinical.rda"))
s.patient <-  c("bcr_patient_barcode", "gender","vital_status","days_to_birth", "days_to_death", "days_to_last_followup","race_list", "stage_event_pathologic_stage")

s.df <- NA
for (i in s.patient) {
    s.df <- cbind(s.df, patient[grep(i, colnames(patient))])
}
s.df <- s.df[,-1]
df_patient <- s.df[!duplicated(s.df),]

# Process patient info
df_patient$vital_status <- sapply(as.numeric(df_patient$vital_status), function(x) x - 1)
df_patient <- df_patient %>% mutate(survival_time = coalesce(days_to_death, days_to_last_followup)) # there are some negative days
max.censored <- max(df_patient$survival_time[df_patient$vital_status == 0])
df_patient <- df_patient[df_patient$survival_time >= 0,]
df_patient$vital_status[df_patient$survival_time == max.censored] <- 1
df_patient$imputed.log.times <- impute.survival(df_patient$survival_time, df_patient$vital_status)
df_patient$age <- floor(df_patient$days_to_birth/365)/-1# Convert age

extra <- c("days_to_birth", "survival_time", "days_to_last_followup", "days_to_death")
df_patient <- dplyr::select(df_patient, -extra)
df_patient <- df_patient[df_patient$imputed.log.times >= 0,]

df_patient[df_patient==""] <- NA

for (c in colnames(df_patient)) {

    if (!is.numeric(df_patient[,c]) && c != "bcr_patient_barcode") {
        which.one <- which( levels(df_patient[,c]) == "")
        levels(df_patient[,c])[which.one] <- NA
        df_patient[,c] = sapply(sapply(df_patient[,c], as.factor), as.numeric) 
        print(paste0(c, " is altered")) 
    }
}

df_patient = na.omit(df_patient)
print("Processed patient dataframe: ")
head(df_patient)

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
                    save.filename = paste0(project, "_exp.rda"))
prep <- GDCprepare(g_query) 
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
exp_matrix$patient = NULL
rownames(exp_matrix) <- exp_matrix$bcr
exp_matrix <- dplyr::select(exp_matrix, -c(bcr, patient))

for (c in c(1:4)) {

which.one <- which( levels(exp_matrix[,c]) == "")
levels(exp_matrix[,c])[which.one] <- NA
print(paste0(c, "is altered"))
exp_matrix[,c] = sapply(sapply(exp_matrix[,c], as.factor), as.numeric) 
}


#Format patient
tmp = strsplit(rownames(exp_matrix), "-")
tmp = unlist(lapply(tmp, function(x) paste(x[[1]], x[[2]], x[[3]], sep = "-")))
rownames(exp_matrix) <- unlist(tmp)

exp_matrix$donorId <- rownames(exp_matrix)
exp_matrix = dplyr::select(exp_matrix, -bcr)

# Select for DEEG only
cancer_DEG = read.csv(paste0("../DEA/", project, "_DEG.csv"))
cancer_DEG = cancer_DEG$x
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
DEG_ensmbl = getBM(attributes = "ensembl_gene_id", filters = "hgnc_symbol", values = cancer_DEG, mart = ensembl)
DEG_ensmbl = DEG_ensmbl$ensembl_gene_id
intersect_DEG = DEG_ensmbl[DEG_ensmbl %in% colnames(exp_matrix)]
exp_matrix = dplyr::select(exp_matrix, c("donorId", "TSS", "portion", "plate", "center", intersect_DEG))


# 4. Prepare covariate matrix, whole dataset matrix and a vectoor of treatment types
colnames(df_patient)[1] = "donorId"
whole_dataset = inner_join(df_patient, as.data.frame(exp_matrix), by = "donorId")
whole_dataset = dplyr::select(whole_dataset, -vital_status)

# This renaming step is critical as the HTE main function will rely on the column named outcome to indicate Y
colnames(whole_dataset)[colnames(whole_dataset) == "imputed.log.times"] = "outcome"
covar_mat= dplyr::select(whole_dataset, -c("donorId", "outcome"))
tx_vector = intersect_DEG

obsNumber <- dim(covar_mat)[1]
trainId <- sample(1: obsNumber, floor(obsNumber/2), replace = FALSE)
registerDoParallel(10)

result <- run.hte(covar_mat, tx_vector, whole_dataset, project, covar_type = "expression", trainId, seed = 111, is.binary = T, is_save = T, save_split = T, is.tuned = F, thres = 0.75, n_core = 8, output_directory = output_file)
write.csv(result[[1]], paste0(output_file, '_expression_correlation_test_result.csv'), quote = F, row.names = F)
write.csv(result[[2]], paste0(output_file, '_expression_calibration_result.csv'), quote = F, row.names = F)
write.csv(result[[3]], paste0(output_file, '_expression_median_t_test_result.csv'), quote = F, row.names = F)
write.csv(result[[4]], paste0(output_file, '_expression_permutate_testing_result.csv'), quote = F, row.names = F)