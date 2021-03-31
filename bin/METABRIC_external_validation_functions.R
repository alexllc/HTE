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
#' @param cf forest object from `grf::causal_forest`
#' @param newCovar the external matrix from an OOB dataset with the same number of columns as the covariate matrix used to build the `cf`. Default is `NULL`, indicating an OOB test is not required.
#' @param save (saving option slot reserved for future update)
#' @param OriginalTauPred the in-bag tau values generated for the patient entries in `newCovar`, used for comparison with the tau prediction of an OOB forest. Default is NULL, indicating an OOB test is not required. Does not have an effect if `OriginalTauPred` is not `NULL` but `OriginalTauPred` is `NULL`
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
        lm_test <- lm(pred[,1] ~ OriginalTauPred + 0)
        print(lm_test)
        return(c(overall_tau_simes_p, test_corr, mean_test_res$statistic, mean_test_res$p.value, mean_test_res$conf.int, lm_test$coefficients))
    }
    return(c(overall_tau_simes_p, tc[1,1], tc[1,4], tc[2,1], tc[2,4]))
}

#' Function for retreiving CNA table from METABRIC dataset
#' 
#' @param
#' 
fetch_cna <- function() {
    if( !file.exists("./raw/cbioportal_METABRIC/data_CNA.txt") )
        message("Please download METABRIC file and exttact to ./raw/cbioportal_METABRIC")
    
    cna <- fread("./raw/cbioportal_METABRIC/data_CNA.txt")
    cna = as.data.frame(t(as.matrix(cna)))
    colnames(cna) = cna[1,]
    cna = cna[-c(1,2),]
    cna = cna[complete.cases(cna),]
    cna = sapply(cna, as.numeric)
    cna = apply(cna, 2, sum)
    return(cna)
}

#' Function to print out summary data of two CF
#' 
#' @param cf1 forest object from `grf::causal_forest`
#' @param cf2 forest object from `grf::causal_forest`
#' @param covar1 the covar matrix used to train cf1
#' @param covar2 the covar matrix used to train cf2
#' @param fit_dat_index the covar used to fit (1 or 2)
#' @param save (saving option slot reserved for future update)
#' 
corr_test_cf <- function(cf1 = NULL,
                         cf2 = NULL,
                         covar1 = NULL,
                         covar2 = NULL,
                         fit_dat_index = 1,
                         save = FALSE) {

    if (fit_dat_index == 1) {
        # Make prediction
        pred_cf1_covar1 <- predict(cf1, newdata = covar1, estimate.variance = TRUE, num.threads = 8)$predictions
        pred_cf2_covar1 <- predict(cf2, newdata = covar1, estimate.variance = TRUE, num.threads = 8)$predictions

        # # Perform best_linear_projection() to check whether a same dataset can show similar conditional treatment effect
        # blp_cf1 <- best_linear_projection(cf1, A = covar1)
        # blp_cf2 <- best_linear_projection(cf2, A = covar1)

        # Perform data consistency test if new covarite matrix from another dataset is given
        # Aim at testing whether two model will give similar prediction on the same dataset, so called external validation.
        ke_test <- ks.test(pred_cf1_covar1, pred_cf2_covar1) # two-sample Kolmogorov-Smirnov test to see whether prediction are from same distribution
        test_corr <- correlation_test(pred_cf1_covar1, pred_cf2_covar1, methods = c("pearson", "kendall", "spearman"), alt = "two.sided")
        mean_test_res <- SIGN.test(pred_cf1_covar1, pred_cf2_covar1)

        lm_test <- summary(lm(pred_cf1_covar1 ~ pred_cf2_covar1 + 0))
        # Calculate 95% CI
        lmlci <- lm_test$coefficients[1] - 1.96 * lm_test$coefficients[2]
        lmuci <- lm_test$coefficients[1] + 1.96 * lm_test$coefficients[2]

        deming_test <- mcr::mcreg(x = pred_cf2_covar1, y = pred_cf1_covar1, method.reg = "Deming", method.ci = "jackknife")@para

        # return vector
        # Kolmogorov-Smirnov test pval; correlation_test; mean_test_res statistic; mean_test_res p value; mean_test_res confident interval; linear regression coeff; linear regression coeff se; linear regression coeff pval; deming reg slope; deming reg slope se; deming reg slope 95% CI
        return(c(ke_test$p.value, test_corr, mean_test_res$statistic,
                mean_test_res$p.value, mean_test_res$conf.int,
                lm_test$coefficients[1], lm_test$coefficients[2], paste0("[", lmlci, ", ", lmuci, "]"), lm_test$coefficients[4], deming_test[2], deming_test[4], paste0("[", deming_test[6], ", ", deming_test[8], "]")))
    } else {
        # Make prediction
        pred_cf1_covar2 <- predict(cf1, newdata = covar2, estimate.variance = TRUE, num.threads = 8)$predictions
        pred_cf2_covar2 <- predict(cf2, newdata = covar2, estimate.variance = TRUE, num.threads = 8)$predictions

        # # Perform best_linear_projection() to check whether a same dataset can show similar conditional treatment effect
        # blp_cf1 <- best_linear_projection(cf1, A = covar2)
        # blp_cf2 <- best_linear_projection(cf2, A = covar2)

        # Perform data consistency test if new covarite matrix from another dataset is given
        # Aim at testing whether two model will give similar prediction on the same dataset, so called external validation.
        ke_test <- ks.test(pred_cf1_covar2, pred_cf2_covar2) # two-sample Kolmogorov-Smirnov test to see whether prediction are from same distribution
        test_corr <- correlation_test(pred_cf1_covar2, pred_cf2_covar2, methods = c("pearson", "kendall", "spearman"), alt = "two.sided")
        mean_test_res <- SIGN.test(pred_cf1_covar2, pred_cf2_covar2)

        lm_test <- summary(lm(pred_cf1_covar2 ~ pred_cf2_covar2 + 0))
        # Calculate 95% CI
        lmlci <- lm_test$coefficients[1] - 1.96 * lm_test$coefficients[2]
        lmuci <- lm_test$coefficients[1] + 1.96 * lm_test$coefficients[2]

        deming_test <- mcr::mcreg(x = pred_cf2_covar2, y = pred_cf1_covar2, method.reg = "Deming", method.ci = "jackknife")@para

        # return vector
        # Kolmogorov-Smirnov test pval; correlation_test; mean_test_res statistic; mean_test_res p value; mean_test_res confident interval; linear regression coeff; linear regression coeff se; linear regression coeff pval; deming reg slope; deming reg slope se; deming reg slope 95% CI
        return(c(ke_test$p.value, test_corr, mean_test_res$statistic,
                mean_test_res$p.value, mean_test_res$conf.int,
                lm_test$coefficients[1], lm_test$coefficients[2], paste0("[", lmlci, ", ", lmuci, "]"), lm_test$coefficients[4], deming_test[2], deming_test[4], paste0("[", deming_test[6], ", ", deming_test[8], "]")))
    }

}

#' Function for performing best linear projection
#' 
#' @param cf the causal forest model
#' @param covar the covariates we want to project the CATE onto
#' 

best_linear_proj_func <- function(cf = NULL,
                                  covar = NULL) {

    best_proj <- grf::best_linear_projection(forest = cf, A = covar)

    best_proj_df <- as.data.frame(matrix(best_proj, dim(best_proj)[1]))
    colnames(best_proj_df) <- c("estimate", "std_Error", "t_value", "p-value")
    row.names(best_proj_df) <- c("Intercept", colnames(covar))
    best_proj_df <- best_proj_df[best_proj_df$`p-value` < 0.05, ] # Select the variable whose estimate is siginifcant.

    return(best_proj_df)
}