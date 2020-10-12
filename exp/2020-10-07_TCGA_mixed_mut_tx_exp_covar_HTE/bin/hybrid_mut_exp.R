# Script to run an HTE analysis using gene mutation as treatment and expression as covariates

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

# Experimentation
library(varSelRF)
library(VSURF)

# Set cancer type
# (testing on HNSC)
project = "HNSC"


# Load HTE function sources
source("../grf_parameters.R")
source("../HTE_main_functions.R")
source("../HTE_validation_functions.R")
source("../NNMIS_survival_imputation.R")

####################### Prepare clinical information #######################

if (!file.exists("../TCGA_CDR_clean.csv")) {

    download.file("https://ars.els-cdn.com/content/image/1-s2.0-S0092867418302290-mmc1.xlsx", "TCGA-CDR-SupplementalTableS1.xlsx")

    cdr = read_excel("../TCGA-CDR-SupplementalTableS1.xlsx")
    clinical_dat = dplyr::select(cdr, c(bcr_patient_barcode, type, age_at_initial_pathologic_diagnosis,  gender, ajcc_pathologic_tumor_stage, tumor_status, OS, OS.time, tumor_status))

    clinical_dat[clinical_dat == "#N/A"] <- NA
    clinical_dat <- subset(clinical_dat, !is.na(OS.time) & !is.na(OS) & !is.na(age_at_initial_pathologic_diagnosis))
    colnames(clinical_dat)[colnames(clinical_dat)=="bcr_patient_barcode"] = "donorId"
    clinical_dat = as.data.frame(clinical_dat)
    write.csv(clinical_dat, "../TCGA_CDR_clean.csv", row.names = F)
} else {
    clinical_dat = read.csv("../TCGA_CDR_clean.csv")
}


# specify cancer type here
ss_patient <- dplyr::filter(clinical_dat, type == project)
ss_patient = impute_with_NNMIS(ss_patient)
ss_patient = dplyr::select(ss_patient, -c(OS, OS.time))
ss_patient = ss_patient[complete.cases(ss_patient),]
print("Processed patient dataframe: ")

head(ss_patient)

####################### Load mutation dataset #######################
mut = read.csv(paste0("../../mut_HTE/", project, "_all_mut_freq.csv"))

####################### Load expression covariates #######################
if (!file.exists(paste0("../../expression_HTE/HTSeqData/", project, "_exp.rda")) ) {
    
    # Fetch expression data from GDC
    g_query <- GDCquery(project = paste0("TCGA-", project),
                    data.category = "Transcriptome Profiling",
                    legacy = F,
                    data.type = "Gene Expression Quantification",
                    workflow.type = "HTSeq - FPKM-UQ",
                    sample.type = "Primary Tumor")

    GDCdownload(g_query)
    expdat <- GDCprepare(query = g_query,
                        save = TRUE,
                        save.filename = paste0("./HTSeqData/",project, "_exp.rda"))
    prep <- GDCprepare(g_query) 
} else {
    load(paste0("../../expression_HTE/HTSeqData/", project, "_exp.rda"))
    prep = data
}

exp_matrix <- SummarizedExperiment::assay(prep, "HTSeq - FPKM-UQ")

# Only select primary tumor samples
# Somehow this line removed so many patients
exp_matrix <- exp_matrix[,grep("-01.-", colnames(exp_matrix))]

#Scale FPKM value by patients
exp_matrix <- exp_matrix/apply(exp_matrix,2,max)
exp_matrix <- t(exp_matrix)

# Address batch effects by extending tissue source, aliquot, plate and sequencing center as additional covariates
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


####################### Running an HTE for each gene, using the same expression dataset as covariate matrix #######################

# Make sure all matrcies do not have NA values or GRF will not accept this input
exp_matrix = exp_matrix[complete.cases(exp_matrix),]
clinical_dat = clinical_dat[!is.na(clinical_dat$tumor_status),]


# Select common patients between expresesion, mutation and clinical datasets
common_pat = intersect(exp_matrix$donorId, clinical_dat$donorId)
common_pat = intersect(common_pat,  mut$donorId)

# Map index of common_pat with clinical_data and select only the common patient rows
pat_index = sapply(common_pat, function(x) which(clinical_dat$donorId == x))
Y = clinical_dat$tumor_status[pat_index] 
Y = as.numeric(Y == "WITH TUMOR") # convert outcome event into numeric
pat_index = sapply(common_pat, function(x) which(exp_matrix$donorId == x))
X = exp_matrix[pat_index,]
X$donorId = NULL # located at the end of the matrix


# Set treatment genes
tx_gene_ls = "TP53" # we test this approach by the popular genes first


for (tx_gene in tx_gene_ls) {
    pat_index = sapply(common_pat, function(x) which(mut$donorId == x))
    tx = mut[pat_index,tx_gene]

    # Build boosted regression forest to pick out the relevant covariates, the commented out parameters are default
    boosted_rf = boosted_regression_forest(
                    X,
                    Y,
                    num.trees = 5000, 
                    seed = 2020# increased to 5k
                    # sample.weights = NULL,
                    # clusters = NULL,
                    # equalize.cluster.weights = FALSE,
                    # sample.fraction = 0.5,
                    # mtry = min(ceiling(sqrt(ncol(X)) + 20), ncol(X)),
                    # min.node.size = 5,
                    # honesty = TRUE,
                    # honesty.fraction = 0.5,
                    # honesty.prune.leaves = TRUE,
                    # alpha = 0.05
                    # imbalance.penalty = 0,
                    # ci.group.size = 2,
                    # tune.parameters = "none",
                    # tune.num.trees = 10,
                    # tune.num.reps = 100,
                    # tune.num.draws = 1000,
                    # boost.steps = NULL,
                    # boost.error.reduction = 0.97,
                    # boost.max.steps = 5,
                    # boost.trees.tune = 10,
                    # num.threads = NULL,
                )
    boosted_predict = predict(boosted_rf)

    # Since we have 50k genes, we will try to pre-select genes by using varSelRF a new package I found on Utah state U lecture
    set.seed(2020)
    # varsel = varSelRF(X, as.factor(Y))
    # Using 
    varsel_ranger = VSURF(X, Y, RFimplem = "ranger", parallel = TRUE, verbose = TRUE)

    X_subset = X[,varsel$selected.vars]
    cf = causal_forest(X_subset, Y, tx)
    cf_predict = predict(cf, estimate.variance = TRUE)
    cf_stats = compute_stats(cf_predict)

    # CF without using the preselection
    cf_gen = causal_forest(X, Y, tx)
    cf_predict_gen = predict(cf_gen, estimate.variance = TRUE)
    cf_stats_gen = compute_stats(cf_predict_gen)
}