#' Script to perform external validation with METABRIC 
# Run from base directory
setwd("~/project/HTE/")

# source bin scripts
bin_ls = list.files("./bin")
for (bin in bin_ls){
    source(paste0("./bin/", bin))
}

# Parameters for this run
validis_tuned = FALSE
is_binary = TRUE
mut_thres = 0.01 # 0.05 yeilded too little results


# Download METABRIC and TCGA files from cBioportal if not already done so, these two conditions prevents script progression if raw files not available, you will need to extract 
if( !file.exists("./raw/cbioportal_METABRIC/")) {
    download.file("http://download.cbioportal.org/brca_metabric.tar.gz", "./raw/")
    untar("./raw/brca_metabric.tar.gz", exdir = "./dat/METABRIC-TCGA_external_validation/")
}

if (!file.exists("./dat/METABRIC-TCGA_external_validation/")){
    download.file("http://download.cbioportal.org/brca_tcga_pan_can_atlas_2018.tar.gz", "./raw/")
    untar("./raw/‘brca_tcga_pan_can_atlas_2018.tar.gz.1’", exdir = "./dat/METABRIC-TCGA_external_validation/")
}

# I. PREPARE FROM RAW DATASETS
# 1. Fetch clinical data
metab_clin <- fetch_metab_clinical() # impute with simple KM method, OS in months
metab_clin <- mk_id_rownames(metab_clin)
tcga_clin <- fetch_clinical_data("BRCA", outParam = "OS", col_vec = c("bcr_patient_barcode", "type", "age_at_initial_pathologic_diagnosis", "OS", "OS.time"), outUnitDays2Month = TRUE) # tyep is still required here to subset appropriate patients from the TCGA-CDR. Also, convert to months to match with metabric data.
tcga_clin$type <-NULL
tcga_clin <-mk_id_rownames(tcga_clin)

# 2. Fetch mutation data, note the colnames corresponding to patient id for both dataset are different
metab_mut <-fetch_metab_mut()
metab_mut <-mk_id_rownames(metab_mut)
tcga_mut <-as.data.frame(fread("./dat/mutation_frequencies/MAF_wide_BRCA.csv.gz")) # previously generated, regenerate with fetch_mut_data if not available
tcga_mut <-mk_id_rownames(tcga_mut)

# 3. Fetch expresion median data
metab_z <-fetch_mrna_z_score("metab", save = TRUE)
tcga_z <-fetch_mrna_z_score("tcga", save = TRUE)
tcga_z <- mk_id_rownames(tcga_z) # first colname V1
rownames(tcga_z) <- format_tcga_patient(rownames(tcga_z)) # z score data set included the sample type (i.e. "-01") even though only one type is present in the BRCA dataset
tcga_z$donorId <- NULL # still needs to be removed
metab_z <-mk_id_rownames(metab_z)

# II. HOMOGENEIZE THE TWO DATASETS
# Subset patients with all: clinical, mutation and expression record and perform subsetting with rownames
metab_common_pat <- intersect(rownames(metab_clin), rownames(metab_mut))
metab_common_pat <- intersect(metab_common_pat, rownames(metab_z))
metab_clin <- metab_clin[metab_common_pat,]
metab_mut <- metab_mut[metab_common_pat,]
metab_z <- metab_z[metab_common_pat,]

tcga_common_pat <- intersect(rownames(tcga_clin), rownames(tcga_mut))
tcga_common_pat <- intersect(tcga_common_pat, rownames(tcga_z))
tcga_clin <- tcga_clin[tcga_common_pat,]
tcga_mut <- tcga_mut[tcga_common_pat,]
tcga_z <- tcga_z[tcga_common_pat,]

# Manually add suffix to both mutation and expression matrices
colnames(metab_mut) <- paste0(colnames(metab_mut), ".m")
colnames(metab_z) <- paste0(colnames(metab_z), ".z")

colnames(tcga_mut) <- paste0(colnames(tcga_mut), ".m")
colnames(tcga_z) <- paste0(colnames(tcga_z), ".z")

# Save list of common expression panel genes in both datasets
common_z <-intersect(colnames(metab_z), colnames(tcga_z))
common_mut <- intersect(colnames(metab_mut), colnames(tcga_mut))

# Count the number of patients who have a mutation in these genes
# TGCA
freq_sum <- sapply(tcga_mut, function(x) sum(x != 0))
mut_selection <- freq_sum > dim(tcga_mut)[1]*mut_thres # select only genes with enough observations
tx_list <- colnames(tcga_mut)[mut_selection]

# METABRIC
freq_sum <- sapply(metab_mut, function(x) sum(x != 0))
mut_selection <- freq_sum > dim(metab_mut)[1]*mut_thres # select only genes with enough 
tx_list <- intersect(tx_list, colnames(metab_mut)[mut_selection])
tx_list <- sapply(strsplit(tx_list,"\\."), function(x) x[[1]])

# Prepare covariate matrices for both datasets
metab_X <- merge(metab_mut, metab_z[,common_z], by = 0, all = FALSE, suffixes = c("",""))
metab_X <- mk_id_rownames(metab_X)
metab_X <- merge(metab_clin, metab_X, by = 0, all = FALSE, suffimetab_Xes = c("", ""))
metab_X <- mk_id_rownames(metab_X)
metab_X$outcome <- NULL
metab_X <- dplyr::select(metab_X, all_of(c("AGE_AT_DIAGNOSIS", common_mut, common_z)))

tcga_X <- merge(tcga_mut, tcga_z[,common_z], by = 0, all = FALSE, suffixes = c("",""))
tcga_X <- mk_id_rownames(tcga_X)
tcga_X <- merge(tcga_clin, tcga_X, by = 0, all = FALSE, suffixes = c("", ""))
tcga_X <- mk_id_rownames(tcga_X)
tcga_X$outcome <- NULL
tcga_X <- dplyr::select(tcga_X, all_of(c("age_at_initial_pathologic_diagnosis", common_mut, common_z)))

# Initialize result dataframes
ext_valid_tcga_cf <- data.frame()
ext_valid_metab_cf <- data.frame()

ext_valid_rn <- c("gene", "original_simes_pval", "mean_forest_pred_coeff", "mean_forest_pred_pval", "diff_forest_pred_coeff", "diff_forest_pred_pval", "oob_tau_simes_pval", "pearson_coeff", "pearson_pval", "kendall_coeff", "kendall_pval", "spearman_coeff", "spearman_pval", "SIGN_statistic", "SIGN_pval", "SIGN_lower_CI", "SIGN_upper_CI", "lm_Y_x_coeff")

for(tx in tx_list) {
    message(paste0("Processing gene: ", tx))
    # Build individual forest

    # Check if expression of that gene is also included in the covariate matrix, if it does, remove from the covaraite matrix as well
    metab_W <- as.numeric(metab_X[,which(colnames(metab_X) == paste0(tx, ".m"))] > 0)
    if( paste0(tx, ".z") %in% colnames(metab_X)) {
        tmp_metab_X <- dplyr::select(metab_X, -c(paste0(tx, ".m"), paste0(tx, ".z")))
    } else {
       tmp_metab_X <- dplyr::select(metab_X, -paste0(tx, ".m"))
    }
    metab_cf <-causal_forest(
                            X = tmp_metab_X,
                            Y = metab_clin$outcome,
                            W = metab_W,
                            num.trees = 5000
                        )
    # make prediction for metab forest
    message("Built METABRIC forest: ")
    original_metab_cf_sum <- print_cf_sum(cf = metab_cf)

    tcga_W <- as.numeric(tcga_X[,which(colnames(tcga_X) == paste0(tx, ".m"))] > 0)
    if ( paste0(tx, ".z") %in% colnames(tcga_X)) {
        tmp_tcga_X <- dplyr::select(tcga_X, -c(paste0(tx, ".m"), paste0(tx, ".z")))
    } else {
       tmp_tcga_X <- dplyr::select(tcga_X, -paste0(tx, ".m"))
    }
    tcga_cf <-causal_forest(
                            X = tmp_tcga_X,
                            Y = tcga_clin$outcome,
                            W = tcga_W,
                            num.trees = 5000
                        )
    message("Built TCGA forest: ")
    # make prediction for tcga forest
    original_tcga_cf_sum <- print_cf_sum(cf = tcga_cf)

    message("Cross over OOB validation: ")
    # Corssing over between the two forests
    tcga_in_metab <- print_cf_sum(cf = metab_cf, X = tmp_metab_X, tx = tx, newCovar = tmp_tcga_X, OriginalTauPred = predict(tcga_cf)[,1]) # fit tcga covar into metab cf
    ext_valid_tcga_cf <- rbind(ext_valid_tcga_cf, c(original_metab_cf_sum, tcga_in_metab)) # append to summary result table
    
    metab_in_tcga <- print_cf_sum(cf = tcga_cf, X = tmp_tcga_X, tx = tx, newCovar = tmp_metab_X, OriginalTauPred = predict(metab_cf)[,1]) # fit metab covar into tcga cf
    ext_valid_metab_cf <- rbind(ext_valid_metab_cf, c(original_tcga_cf_sum, metab_in_tcga))
}

ext_valid_tcga_cf <- cbind(tx_list, ext_valid_tcga_cf)
colnames(ext_valid_tcga_cf) <- ext_valid_rn
# calcualte calibration lm z score for tcga
ext_valid_tcga_cf$lm_Z <- abs(ext_valid_tcga_cf$lm_Y_x_coeff - 1)/sd(ext_valid_tcga_cf$lm_Y_x_coeff)
write.csv(ext_valid_tcga_cf, "./exp/2020-10-27-external_validation_mixed_HTE/res/TCGA_METABRIC-BRICA_external_validtaion_TCGA_forest.csv", row.names = FALSE)

# calcualte calibration lm z score for metab
ext_valid_metab_cf <- cbind(tx_list, ext_valid_metab_cf)
colnames(ext_valid_metab_cf) <- ext_valid_rn
ext_valid_metab_cf$lm_Z <- abs(ext_valid_metab_cf$lm_Y_x_coeff - 1)/sd(ext_valid_metab_cf$lm_Y_x_coeff)
write.csv(ext_valid_METAB_cf, "./exp/2020-10-27-external_validation_mixed_HTE/res/TCGA_METABRIC-BRICA_external_validtaion_METABRIC_forest.csv", row.names = FALSE)