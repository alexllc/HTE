# Script for running HTE analysis using IQLR transformed gene expression as both treatment and covaraites

## Run from base directory
setwd("~/project/HTE/")

## source bin scripts
bin_ls = list.files("./bin")
for (bin in bin_ls){
    source(paste0("./bin/", bin))
}

## Set parameters for this run
cancer_type <- "BRCA"
endpt <- "OS"
output_file <- "./exp/2021-01-15_ALDEx_DEA/res/"
resume <- 1115
useALDEX2_DE <- TRUE
use_DE_covar_only <- TRUE

## Prepare clinical dataframe
clinical <- fetch_clinical_data(cancer_type, outParam = endpt, imputeMethod = "simple", outUnitDays2Month = TRUE, discard = c("type", "tumor_status"))
clinical <- mk_id_rownames(clinical)

## Load expression matrix
if(!DE_covar_only) {
    exp <- as.data.frame(fread("./dat/BRCA_iqlr_expected_count.csv.gz"))
    rownames(exp) <- exp$V1
    exp$V1 <- NULL
    exp <- t(exp)
} else {
   if (!useALDEX2) {

   }
}
# select primary tumor only
CancerProject <- "TCGA-BRCA"

query <- GDCquery(project = CancerProject,
                data.category = "Transcriptome Profiling",
                data.type = "Gene Expression Quantification", 
                workflow.type = "HTSeq - Counts")

samplesDown <- getResults(query,cols=c("cases"))

dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                typesample = "TP")

rownames(exp) <- gsub("\\.", "-", rownames(exp))
# filter replicates
exp <- exp[dataSmTP, ]
selection <- filter_replicate_samples(rownames(exp))
exp <- exp[selection,]

# format patient ID
rownames(exp) <- format_tcga_patient(rownames(exp))

## Subset patients and settings for HTE
# Check the list of common patients across three data frames
common_pat = rownames(clinical)[rownames(clinical) %in% rownames(exp)]
message(paste0("Number of patients with expression and clinical entries is: ", length(common_pat)))

# Outcome vector for causal forest
Y = clinical[common_pat, "outcome"]
X = exp[common_pat,]

# Load DEA results
if (!useALDEX2_DE) {
    dea <- read.csv("./dat/tables/BRCA_DEG_rerun.csv")
    rownames(dea) <- dea$X
    dea$X <- NULL
    dirct <- dea$logFC > 0
    names(dirct) <- dea$mRNA
} else {
   aldex <- readRDS("./dat/ALDEx2_result.rds.gz")
   sig <- filter(aldex, abs(effect) > 1.5 & abs(overlap) < 0.05 ) # not default filter setting
   head(sig)[, c(1:3, 1219:1226)]
    dea <- sig$effect
    names(dea) <- rownames(sig)
}

tx_list = names(dirct)
X = inner_join(rownames_to_column(clinical[common_pat, ]), rownames_to_column(as.data.frame(X)), by = "rowname")

whole_dat <- X
X <- dplyr::select(X, -all_of(c("rowname", "outcome")))
colnames(whole_dat)[colnames(whole_dat) == "rowname"] = "donorId"

# Assign treatment group
W = create_tx_matrix(txVector = tx_list, 
                    binaryVector = rep(TRUE, length(tx_list)), 
                    cutoffThreshDf = data.frame(
                                        dirct = ifelse(dirct, ">", "<"), 
                                        thresh = ifelse(dirct, 0.75, 0.25)
                                        ), 
                    covarMat = X)

# options for resuming
tx_list <- tx_list[resume:length(tx_list)]
W <- W[,resume:length(tx_list)]


## Run and save HTE
obsNumber <- dim(X)[1]
trainId <- sample(1: obsNumber, floor(obsNumber/2), replace = FALSE)
registerDoParallel(10)

result <- run.hte(covar_mat = X, 
                tx_vector = tx_list, 
                whole_dataset = whole_dat, 
                project = cancer_type, 
                W_matrix = W, 
                trainId = trainId, 
                seed = 111, 
                is_binary = T, 
                is_save = T, 
                save_split = T, 
                is_tuned = F, 
                thres = 0.75, 
                n_core = 70, 
                output_directory = output_file, 
                perm_all = FALSE) # pre-filtered thres, should not have an effect here

write.csv(result[[1]], paste0(output_file, cancer_type, '_mixed_HTE_correlation_test_result.csv'), quote = F, row.names = F)
write.csv(result[[2]], paste0(output_file, cancer_type, '_iqlr_calibration_result.csv'), quote = F, row.names = F)
write.csv(result[[3]], paste0(output_file, cancer_type, '_iqlr_median_t_test_result.csv'), quote = F, row.names = F)
write.csv(result[[4]], paste0(output_file, cancer_type, '_iqlr_permutate_testing_result.csv'), quote = F, row.names = F)


# for continuous treatment vectors
cont_result <- run.hte(covar_mat = X, 
                tx_vector = tx_list, 
                whole_dataset = whole_dat, 
                project = cancer_type, 
                W_matrix = X[,tx_list], 
                trainId = trainId, 
                seed = 111, 
                is_binary = F, 
                is_save = T, 
                save_split = T, 
                is_tuned = F, 
                thres = 0.75, 
                n_core = 70, 
                output_directory = "./exp/2021-01-15_ALDEx_DEA/tmp/tmp_res/", 
                perm_all = FALSE) 