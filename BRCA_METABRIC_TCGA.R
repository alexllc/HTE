library(cgdsr)
setwd("../HTE")
project = "BRCA"
source("./Overlap_of_varImp_ver1_asFunction-original.R")
source("./HTE_mutation.R")

aggr_res <- function(res_mat, i, est_col_list = NULL) {
    removed_na <- na.omit(res_mat[, i])
    # simes_pval <- ifelse(length(removed_na) > 0, simes.test(removed_na), NA)
    if(!(i %in% est_col_list)){
        pval <- ifelse(length(removed_na) > 0, simes.partial(2, removed_na), NA)
    } else {
        pval <- ifelse(length(removed_na) > 0, extract_binom_pval(removed_na), NA) 
    }
}

## Import and downloading METABRIC mutation data
mut = fread("./METABRIC/data_mutations_extended.txt")
cnt = mut %>% group_by(Tumor_Sample_Barcode) %>% count(Hugo_Symbol)
cnt = spread(cnt, Hugo_Symbol, n)
cnt[is.na(cnt)] = 0
sel_genes = intersect(colnames(cnt), tx_vector) # should be 167/174 of METABRIC's panel
cnt = dplyr::select(cnt, all_of(c("Tumor_Sample_Barcode", sel_genes)))

## Procress METABRIC clinical data
clinical = fread("./METABRIC/data_clinical_patient.txt", skip = 4, header = T)
clinical = dplyr::select(clinical, PATIENT_ID, OS_MONTHS, OS_STATUS, AGE_AT_DIAGNOSIS)
clinical = na.omit(clinical)
# Convert event to 0/1
clinical$OS_STATUS = as.integer(clinical$OS_STATUS=="1:DECEASED") # METABRIC updated the OS status term
# Convert OS mo to days
clinical$OS_MONTHS = clinical$OS_MONTHS * 30.4167
# Impute survival time
max_censored = max(clinical$OS_MONTHS[clinical$OS_STATUS == 0], na.rm = T)
clinical$OS_STATUS[clinical$OS_MONTHS==max_censored] = 1
clinical$outcome = impute.survival(clinical$OS_MONTHS,  clinical$OS_STATUS)
colnames(clinical)[1] = "Tumor_Sample_Barcode"
clinical$OS_MONTHS = NULL
clinical$OS_STATUS = NULL
head(clinical)

metab_all = inner_join(clinical, cnt, by = "Tumor_Sample_Barcode")
colnames(metab_all)[1] = "donorId"
metab_covar = dplyr::select(metab_all, -c(donorId, outcome))


## Process TCGA data to match METABRIC data
whole_dataset = left_join(ss_patient, dplyr::select(wtcga, all_of(c("donorId", sel_genes))), by = "donorId")
whole_dataset = dplyr::select(whole_dataset, -c("gender", "ajcc_pathologic_tumor_stage"))
whole_dataset = whole_dataset[complete.cases(whole_dataset),]
covar_mat= dplyr::select(whole_dataset, -c("donorId", "outcome"))


## Running HTE for both datasets
covar_type = "mutation"
seed = 111
is.binary = T
is_save = T
save_split = T
is.tuned = F
thres = 0.75
n_core = 8
output_directory = output_file


cf.estimator <- ifelse(is.tuned, cf.tuned, cf)

col_names <- c('simes.pval', 'partial.simes.pval', 'pearson.estimate','pearson.pvalue', 'kendall.estimate','kendall.pvalue', 'spearman.estimate','spearman.pvalue', 'fisher.pval', 't.test.a.pval', 't.test.b.pval')

correlation_test_ret = NULL
corr_test_names = c("gene", "simes.pval", "partial.simes.pval", "pearson.estimate", "pearson.pvalue", "kendall.estimate", "kendall.pvalue", "spearman.estimate", "spearman.pvalue")
correlation_test_ret = NULL 


overlap_test_res = NULL
overlap_test_names = c("gene","fisher_top1","fisher_top3","fisher_top5","fisher_top10","fisher_top20","fisher_top30","above0_est","above0_pval","perm_fisher","perm_min","perm_simes","perm_softomni")

## Skipping TAF1 because tau prediction from this gene is constant
sel_genes = sel_genes[-which(sel_genes == "TAF1")] # 166 for BRCA


for (tx in sel_genes){
    
    print(paste0(c('#', rep('-', 40), ' running ', which(sel_genes == tx), ' of ', length(sel_genes), rep('-', 40)), collapse = ''))
    # tx vector and covar mat for TCGA
    cat('Treatment name:', tx, fill = T)

    T_treatment <- as.data.frame(covar_mat)[, tx]
    T_covariates <- as.matrix(dplyr::select(covar_mat, -tx))
    T_treatment <- as.numeric(T_treatment != 0) # only for mutation
    T_Y <- whole_dataset$outcome

    # tx vector and covar mat for METAB
    M_treatment <- as.data.frame(metab_covar)[, tx]
    M_covariates <- as.matrix(dplyr::select(metab_covar, -tx))
    M_treatment <- as.numeric(M_treatment != 0) # only for mutation
    M_Y <- metab_all$outcome

    if (length(unique(T_treatment) <= 1 | unique(M_treatment)) <= 1) {
        print("Gene mutation distribution too sparse, skipping.")
        next
    }
    file_prefix = paste0(output_directory, project, "_", tx)

    ## Split half
    T_obs <- dim(T_covariates)[1]
    M_obs <- dim(M_covariates)[1]

    correlation_matrix = NULL
    overlap_matrix = NULL

    no_repeats = 10
    for (i in seq(no_repeats)) {
        observation_result_tcga <- matrix(0, nrow = T_obs, ncol = 4)
        observation_result_metab <- matrix(0, nrow = M_obs <- dim(T_covariates)[1], ncol = 4)

        tcga_forest <- cf.estimator(T_covariates, T_Y, T_treatment, seed = seed)
        tau_pred_tcga_train <- predict(tcga_forest, estimate.variance = T, num.threads = n_core)
        tau_pred_tcga_metab <- predict(tcga_forest, newdata = M_covariates, estimate.variance = T, num.threads = n_core)

        metab_forest <- cf.estimator(M_covariates, M_Y, M_treatment, seed = seed)
        tau_pred_metab_train <- predict(metab_forest, estimate.variance = T, num.threads = n_core)
        tau_pred_metab_tcga <- predict(metab_forest, newdata = T_covariates, estimate.variance = T, num.threads = n_core)

        # compute z-score, pvalues, and ajusted.pvalues
        tau_tcga_train_stats <- compute_stats(tau_pred_tcga_train)
        tau_tcga_stats <- compute_stats(tau_pred_tcga_metab)
        tau_metab_train_stats <- compute_stats(tau_pred_metab_train)
        tau_metab_stats <- compute_stats(tau_pred_metab_tcga)

        if(save_split){
            write.csv(tau_tcga_train_stats, file = paste0(file_prefix, '_observation_', i, '_result_tcga.csv'), row.names = F)
            write.csv(tau_metab_stats, file = paste0(file_prefix, '_observation_', i, '_result_tcga_in_metab.csv'), row.names = F)
            
            write.csv(tau_metab_train_stats, file = paste0(file_prefix,  '_observation_', i,'_result_metab.csv'), row.names = F)
            write.csv(tau_tcga_stats, file = paste0(file_prefix, '_observation_', i,'_result_metab_in_tcga.csv'), row.names = F)

            # Extract varimp from each of the split half forest (Jun 13, 2020 @alex)
            T_varImp_extract = variable_importance(tcga_forest, max.depth = 4)
            T_varImp <- data.frame(variable = colnames(T_covariates), T_varImp_extract)
            write.csv(T_varImp, file = paste0(file_prefix, '_observation_', i, '_varimp_tcga.csv'), row.names = F, quote = F)

            M_varImp_extract = variable_importance(metab_forest, max.depth = 4)
            M_varImp <- data.frame(variable = colnames(M_covariates), M_varImp_extract)
            write.csv(M_varImp, file = paste0(file_prefix, '_observation_', i, '_varimp_metab.csv'), row.names = F, quote = F)
            #message(paste0("Varimp for observation ", i, " saved."))
        }

        # Overlap of varimp using Prof So's script (Jun 17, 2020)
        input_matrix = merge(T_varImp, M_varImp, by = "variable")

        overlap = test_overlap_VarImp (input_matrix,
                            no_perm = 5000,
                            no_cluster = 10, #no. of clusters for parallel running 
                            top_percentile_list = c(0.01, 0.03, 0.05, 0.1, 0.2, 0.3) 		
                            )
        overlap_ext = c()
        for (i in 1:length(overlap)) {
            if (i != 2) overlap_ext = c(overlap_ext, overlap[[i]])
        }

        overlap_matrix = rbind(overlap_matrix, overlap_ext)

        simes_pval_tcga <- simes.test(tau_tcga_stats[, 3])
        simes_pval_metab <- simes.test(tau_metab_stats[, 3])

        partial_simes_pval_tcga <- simes.partial(floor(dim(T_covariates)[1] * 0.05), tau_tcga_stats[, 3])
        partial_simes_pval_metab <- simes.partial(floor(dim(M_covariates)[1] * 0.05), tau_metab_stats[, 3])

        # check the correlation between two predictions from two datasets
        test_res_tcga <- correlation_test(tau_tcga_train_stats[, 1], tau_metab_stats[, 1], methods = c('pearson', 'kendall', 'spearman'))
        test_res_metab <- correlation_test(tau_tcga_stats[, 1], tau_metab_train_stats[, 1], methods = c('pearson', 'kendall', 'spearman'))

        fisher_pval_tcga <- fisher.exact.test(tau_tcga_train_stats[, 3], tau_metab_stats[, 3])
        fisher_pval_metab <- fisher.exact.test(tau_tcga_stats[, 3], tau_metab_train_stats[, 3])

        t_test_pval_tcga <- quantile.t.test(tau_tcga_train_stats[, 1], tau_metab_stats[, 1])
        t_test_pval_metab <- quantile.t.test(tau_tcga_stats[, 1], tau_metab_train_stats[, 1]) 
        
        correlation_rslt <- rbind(c(simes_pval_tcga, partial_simes_pval_tcga, test_res_tcga, fisher_pval_tcga, t_test_pval_tcga), c(simes_pval_metab, simes_pval_metab, test_res_metab, fisher_pval_metab, t_test_pval_metab))
        correlation_matrix = rbind(correlation_matrix, correlation_rslt)
    }

    if(is_save){
        colnames(correlation_matrix) <- col_names 
        write.csv(correlation_matrix, file = paste0(file_prefix, '_split_half.csv'), row.names = F, quote = F)
        write.csv(overlap_matrix, file = paste0(file_prefix, '_varimp_overlap.csv'), row.names = F, quote = F)
    }
    aggregated_corr_rslt <- sapply(seq(dim(correlation_matrix)[2]), aggr_res, est_col_list = c(3, 5, 7), res_mat = correlation_matrix)

    aggregated_overlap_rslt <- sapply(seq(dim(overlap_matrix)[2]), aggr_res, est_col_list = 7, res_mat = overlap_matrix)

    ## Alex: Moved to a declared function Jun 27, 2020
    # change to partial_simes_pval 20190929
    # aggregated_rslt <- sapply(seq(dim(correlation_matrix)[2]), function(i){
    #     removed_na <- na.omit(correlation_matrix[, i])
    #     # simes_pval <- ifelse(length(removed_na) > 0, simes.test(removed_na), NA)
    #     if(!(i %in% pval_col_list)){
    #         pval <- ifelse(length(removed_na) > 0, simes.partial(2, removed_na), NA)
    #     } else {
    #         pval <- ifelse(length(removed_na) > 0, extract_binom_pval(removed_na), NA) 
    #     }
    #     return(pval)
    # })



    cat('Fisher extact test pval in trainset:', aggregated_corr_rslt[9], fill = T)
    cat('pearson correlation pval in trainset:', aggregated_corr_rslt[4], fill = T)
    cat('kedall correlation pval in trainset:', aggregated_corr_rslt[6], fill = T)
    cat('spearman correlation pval in trainset:', aggregated_corr_rslt[8], fill = T)

    message("===========Varimp overlap results: ==============")
    for(k in 2:length(aggregated_overlap_rslt)) {
        cat(paste0(colnames(overlap_test_res)[k], ": "), aggregated_overlap_rslt[k], fill = T)
    }

    # Alex @ Jun28, 2020 do.call then append to preset data.frame does not work anymore, switching to rbindlist
    current_ret <- do.call('c', list(list(tx), as.list(aggregated_corr_rslt[1: 8])))
    current_ret <- rbindlist(list(current_ret))
    correlation_test_ret <- rbind(correlation_test_ret, current_ret)

    overlap_tmp <- do.call('c', list(list(tx), as.list(aggregated_overlap_rslt)))
    overlap_tmp <- rbindlist(list(overlap_tmp))
    overlap_test_res <- rbind(overlap_test_res, overlap_tmp)
                             
}

colnames(correlation_test_ret) <- corr_test_names
write.csv(correlation_test_ret, paste0(output_file, project, '_validation_correlation_test_result.csv'), quote = F, row.names = F)
write.csv(overlap_test_res, paste0(output_file, project, '_validation_varimp_overlap_result.csv'), quote = F, row.names = F)