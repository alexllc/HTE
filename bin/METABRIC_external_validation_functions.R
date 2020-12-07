# Script containing functions for performing external validations of TCGA HTE results with METABRIC data.

#' Function for extratcing basic clinical information and imputing survival times with censoring for HTE analysis
#' 
#' @param imputeMethod [string] impute with simple survfit or with NNMIS
#' @param col_vec [string vector] clinical column that needs to be retreived
fetch_metab_clinical <- function(imputeMethod = "simple", col_vec = c("PATIENT_ID", "AGE_AT_DIAGNOSIS", "OS_MONTHS", "OS_STATUS") ) {
    clin = fread("./raw/cbioportal_METABRIC/data_clinical_patient.txt", skip = 4)
    clin = dplyr::select(clin, all_of(col_vec))
    clin = clin[which( !(clin$OS_MONTHS <= 0 | is.na(clin$OS_MONTHS) | clin$OS_STATUS == "")),]
    clin$OS_STATUS = as.integer(clin$OS_STATUS=="1:DECEASED") # METABRIC updated 
    clin$outcome = exp(impute.survival(clin$OS_MONTHS, clin$OS_STATUS)) # no need to take 30.417 reciprical to convert to months here
    clin = dplyr::select(clin, -c(OS_MONTHS, OS_STATUS))
    return(as.data.frame(clin))
}


#' Function to extract mutation frequencies from data_mutations_extended.txt downloaded 
#' 
#' 
#' @param sel_genes (optional)[string vector] list of genes you want to extract from the dataframe, default is extracting all available.
fetch_metab_mut <- function(sel_genes = NULL){
    mut = fread("raw/cbioportal_METABRIC/data_mutations_extended.txt")
    mut = mut %>% group_by(Tumor_Sample_Barcode) %>% count(Hugo_Symbol)
    mut = spread(mut, Hugo_Symbol, n)
    mut[is.na(mut)] = 0
    if(!is.null(sel_genes)) mut = dplyr::select(mut, all_of(c("Tumor_Sample_Barcode", sel_genes)))
    return(as.data.frame(mut))
}


#' Function for extracting mRNA median expression z score from the txt downloaded from cBioPortal
#' 
#' @param data [string] either "metab" or "tcga"
#' @param [bool] save the transposed data frame
fetch_mrna_z_score <- function(data = NULL, save = FALSE) {
        if(!file.exists(paste0("./dat/METABRIC-TCGA_external_validation/exp_median_z/", data, "_mRNA_med_z_scores_transposed.csv.gz"))) {
            if (data == "metab") {
                mrna_z = fread("./raw/cbioportal_METABRIC/data_mRNA_median_Zscores.txt")
            } else {
                mrna_z = fread("./raw/cbioportal_TCGA-BRCA_pancaner_atlas/data_RNA_Seq_v2_mRNA_median_Zscores.txt")
            }
            message(paste0("Transposing ", data, " expression z-score matrix, this could take a *LONG* while."))
            mrna_z = t(mrna_z)
            colnames(mrna_z) = mrna_z[1,]
            mrna_z = mrna_z[-c(1:2),]
            class(mrna_z) = "numeric"
            mrna_z = as.data.frame(mrna_z)
            mrna_z$donorId = rownames(mrna_z)
            if(save) write.csv(mrna_z, file = gzfile(paste0("./dat/METABRIC-TCGA_external_validation/exp_median_z/", data, "_mRNA_med_z_scores_transposed.csv.gz")), row.names = TRUE)
    } else {
        mrna_z = as.data.frame(fread(paste0("./dat/METABRIC-TCGA_external_validation/exp_median_z/", data, "_mRNA_med_z_scores_transposed.csv.gz")))
    }
    z_not_na = sapply(mrna_z, function(x)all(!is.na(x))) # many genes are fully NA
    mrna_z = mrna_z[,z_not_na]
    return(mrna_z)
}

#' Function to only print out summary data of an individual CF
#' 
#' @param cf
#' @param newCovar
#' @param save
#' 
print_cf_sum <- function(cf = NULL,
                        X = NULL,
                        tx = NULL,
                        newCovar = NULL,
                        OriginalTauPred = NULL, 
                        save = FALSE) {
    if(is.null(newCovar)) {
        pred <- predict(cf, estimate.variance = TRUE)
        tc <- test_calibration(cf)
    } else {
       pred <- predict(cf, newdata = newCovar, estimate.variance = TRUE)
    }

    stats <- compute_stats(pred)
    overall_tau_simes_p <- simes.test(stats$tau.pval)
    message(paste0("simes tau p value is: ", overall_tau_simes_p))

    # Perform data consistency test if new covarite matrix from another dataset is given
    if(!is.null(newCovar)) {
        test_corr <- correlation_test(OriginalTauPred, pred[, 1], methods = c("pearson", "kendall", "spearman"), alt = "two.sided")
        print(test_corr)
        mean_test_res <- SIGN.test(OriginalTauPred, pred[,1])
        print(mean_test_res)
        return(c(overall_tau_simes_p, test_corr, mean_test_res$statistic, mean_test_res$p.value, mean_test_res$conf.int))
    }
    return(c(overall_tau_simes_p, tc[1,1], tc[1,4], tc[2,1], tc[2,4]))
}