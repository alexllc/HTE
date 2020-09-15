# library(TCGAbiolinks)
# library(DT)
# library(SummarizedExperiment)
library(plyr)
library(dplyr)
library(tidyr)

#hte_functions
library(grf)
library(MASS)
library(doMC)
library(data.table)
library(survminer)
#library(mppa)
library(doParallel)
library(methods)

# causal inference models
library(grf)
library(BART)
library(ranger)
library(randomForestSRC)
library(randomForest)

# survival imputation
library(survival)

# fetch immune
library(readxl)


#usrwd = "/exeh_4/alex_lau"
usrwd = "/exeh/exe4/alex_lau"
source("./grf_parameters.R")
source("./HTE_main_functions.R")
source("./HTE_validation_functions.R")
# source("./NNMIS_survival_imputation.R")
source("./survival_imputation.R")

cancerType <- "BRCA"#"LUSC" #"LUAD" #"BRCA" #"GBM"

#**********************************************************************************
#                   clinical and mutation freq data from Kai
#**********************************************************************************

setwd("/exeh_3/kai/data/heterogenous_treatment_effects")

#setwd(paste0(usrwd, "/HTE/wd"))
# load TCGA dataset
dataset <- read.csv('TCGA_survival_CDR.csv', header = TRUE, stringsAsFactors = F)
dataset <- dataset[c(2, 5, 4, 32, 33, 3)]   # select relevant information
dataset$gender <- as.numeric(dataset$gender == 'FEMALE')

# some clearning
dataset <- dataset[dataset$PFI.time!=0,]
dataset[dataset == "#N/A"] <- NA
dataset <- subset(dataset, !is.na(PFI.time) & !is.na(PFI) & !is.na(age_at_initial_pathologic_diagnosis))

# load ICGC dataset
icgcMutation <- read.csv('ICGC/ICGC.covatiates.mutation.csv', header = T, stringsAsFactors = F)
donorInfo <- icgcMutation[c(1:4, 8, 5)]
donorInfo$sex <- as.numeric(donorInfo$sex == 'female')

# change names of two dataset
index <- c('donorId', 'gender', 'age', 'censor', 'PFI.time', 'type')
colnames(dataset) <- index
colnames(donorInfo) <- index

# union two datasets
patientInfo <-rbind(dataset, donorInfo)

# specify cancer type here
patientInfo <- subset(patientInfo, type %in% cancerType)
surv.times <- as.numeric(as.character(patientInfo$PFI.time))
cens <- as.numeric(patientInfo$censor)

# get imputed log survival times
max.censored <- max(surv.times[cens == 0])
cens[surv.times == max.censored] <- 1
imputed.log.times = impute.survival(surv.times, cens)

# attach imputed.log.times to original dataset
dataset <- cbind(patientInfo, imputed.log.times)

#**********************************************************************************
#                                mutation data
#**********************************************************************************

mutation.threshold <- 50

# load and process TCGA dataset
load("TCGA_PanCanAtlas_No_Of_Mutations_NoSilentMut_intersected.rda")
donorId <- substr(intersected.mut_adj.rev$Tumor_Sample_Barcode, start=1, stop=12)
mut_adj.rev <- cbind(donorId, intersected.mut_adj.rev[2: dim(intersected.mut_adj.rev)[2]])

# obtain ICGC mutation data
icgcMut <- icgcMutation[c(1, 9: dim(icgcMutation)[2])]
colnames(icgcMut)[1] <- 'donorId'

# rbind two dataset
mutDataset <- rbind(mut_adj.rev, icgcMut)

# merge dataset with mutitation dataset by "bcr_patient_barcode"
merged.dataset <- inner_join(dataset, mutDataset, by="donorId")

# find mutations with mutation count greater than the defined mutation.threhold
pos <- which(colnames(merged.dataset) == "TSPAN6")
freq.mutations <- extract.freq.mutation(merged.dataset, pos, mutation.threshold)

mutation.mat <- merged.dataset[, pos:ncol(merged.dataset)]

# adding confonding variables
age <- merged.dataset$age
gender <- merged.dataset$gender
# type <- as.numeric(merged.dataset$type == 'LUAD')

# remove stage variable first because of lots of missing value
# stage <- as.factor(merged.dataset$clinical_stage)
mutation.mat <- cbind(age, gender, mutation.mat)
mutation.mat$age <- as.numeric(as.character(mutation.mat$age))

#gene_ls = colnames(mutation.mat[pos:ncol(mutation.mat)])
head(mutation.mat[,1:10])


# Further subsetting the treatment genes with driver gene list
# TSG <- read_excel(paste0(usrwd, "/cancerl/data/davoli_genes.xlsx"), sheet = "Table S3A TSG", skip = 2)
# TSG <- head(TSG[order(TSG$TUSON_q_value_TSG, decreasing= F),], n = 300)
# TSG <- TSG$Gene
# OG <- read_excel(paste0(usrwd, "/cancerl/data/davoli_genes.xlsx"), sheet = "Table S3B OG", skip = 2)
# OG <- head(OG[order(OG$TUSON_q_value_OG, decreasing= F),], n = 250)
# OG <- OG$Gene  
# drivers <- c(TSG, OG)

# incl_drivers = drivers[drivers %in% freq.mutations]
#print(paste0('number of mutation will be studied: ', length(incl_drivers) ))

# Comparing with METABRIC
## Only selecting those genes and clinical data (gender) 
freq_sum = sapply(mutation.mat[,3:dim(mutation.mat)[2]], function(x) sum(x != 0))
selection = freq_sum > dim(mutation.mat)[1]*0.05


sel_genes = names(which(selection == T))
mutation.mat = dplyr::select(mutation.mat, c(age, sel_genes))
merged.dataset = dplyr::select(merged.dataset, c(donorId, censor, age, imputed.log.times, sel_genes))
colnames(merged.dataset)[4] = "outcome"


obsNumber <- dim(mutation.mat)[1]
trainId <- sample(1: obsNumber, floor(obsNumber/2), replace = FALSE)


output_file = "/exeh_4/alex_lau/proj/HTE/wd/mut_HTE/rerun_old_code/"
result <- run.hte(covar_mat = mutation.mat, tx_vector = sel_genes, whole_dataset = merged.dataset, project = "BRCA", covar_type = "mutation", trainId = trainId, seed = 111, is.binary = T, is_save = T, save_split = T, is.tuned = F, thres = 0.75, n_core = 8, output_directory = output_file, skip_perm = FALSE)
write.csv(result[[1]], paste0(output_file, cancerType, '_expression_correlation_test_result.csv'), quote = F, row.names = F)
write.csv(result[[2]], paste0(output_file, cancerType, '_expression_calibration_result.csv'), quote = F, row.names = F)
write.csv(result[[3]], paste0(output_file, cancerType, '_expression_median_t_test_result.csv'), quote = F, row.names = F)
write.csv(result[[4]], paste0(output_file, cancerType, '_expression_permutate_testing_result.csv'), quote = F, row.names = F)