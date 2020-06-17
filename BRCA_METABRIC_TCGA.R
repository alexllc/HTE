library(cgdsr)
#setwd("../HTE")

project = "BRCA"
source("./HTE_mutation.R")

## Import and downloading METABRIC data
mut = fread("./METABRIC/data_mutations_extended.txt")
cnt = mut %>% group_by(Tumor_Sample_Barcode) %>% count(Hugo_Symbol)
cnt = spread(cnt, Hugo_Symbol, n)
cnt[is.na(cnt)] = 0
sel_genes = intersect(colnames(cnt), tx_vector) # should be 167/174 of METABRIC's panel
cnt = dplyr::select(cnt, all_of(c("Tumor_Sample_Barcode", sel_genes)))
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
head(clinical)

metab_all = inner_join(clinical, cnt, by = "Tumor_Sample_Barcode")
colnames(metab_all)[1] = "donorId"
metab_covar = dplyr::select(metab_all, -c(donorId, OS_STATUS, outcome))


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

        col_names <- c('simes.pval', 'partial.simes.pval', 'pearson.estimate','pearson.pvalue',
               'kendall.estimate','kendall.pvalue', 'spearman.estimate','spearman.pvalue',
               'fisher.pval', 't.test.a.pval', 't.test.b.pval')

for (tx in sel_genes){

    # tx vector and covar mat for TCGA

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
    for (i in seq(no_repeats)) {
        observation_result.tcga <- matrix(0, nrow = T_obs, ncol = 4)
        observation_result.metab <- matrix(0, nrow = M_obs <- dim(T_covariates)[1], ncol = 4)

        lot = sample(M_obs, T_obs)
        tcga_forest <- cf.estimator(T_covariates, T_Y, T_treatment, seed = seed)
        tau_pred_tcga_train <- predict(tcga_forest, estimate.variance = T, num.threads = n_core)
        tau_pred_tcga_metab <- predict(tcga_forest, newdata = M_covariates, estimate.variance = T, num.threads = n_core)

        metab_forest <- cf.estimator(M_covariates, M_Y, M_treatment, seed = seed)
        tau_pred_metab_train <- predict(metab_forest, estimate.variance = T, num.threads = n_core)
        tau_pred_metab_tcga <- predict(tau.forest.test, newdata = T_covariates, estimate.variance = T, num.threads = n_core)

        # compute z-score, pvalues, and ajusted.pvalues
        tau.a.train.stats <- compute_stats(tau.pred.train.a)
        tau.a.stats <- compute_stats(tau.pred.test.a)
        tau.b.train.stats <- compute_stats(tau.pred.train.b)
        tau.b.stats <- compute_stats(tau.pred.test.b)

        if(save_split){
            observation_result.a[trainId,] <- tau.a.train.stats
            observation_result.a[-trainId,] <- tau.a.stats
            write.csv(observation_result.a, file = paste0(file_prefix,   '_observation_', i, '_result_a.csv'))
            observation_result.b[-trainId,] <- tau.b.train.stats 
            observation_result.b[trainId,] <- tau.b.stats
            write.csv(observation_result.b, file = paste0(file_prefix,  '_observation_', i,'_result_b.csv'))
            
            # Extract varimp from each of the split half forest (Jun 13, 2020 @alex)
            varImp = variable_importance(tau.forest.train)
            train.varimp <- data.frame(variable = varimp_names, varImp, max.depth = 4)
            write.csv(train.varimp, file = paste0(file_prefix, '_observation_', i, '_varimp_train.csv'), row.names = F, quote = F)

            varImp = variable_importance(tau.forest.test)
            test.varimp <- data.frame(variable = varimp_names, varImp, max.depth = 4)
            write.csv(test.varimp, file = paste0(file_prefix, '_observation_', i, '_varimp_test.csv'), row.names = F, quote = F)
            #message(paste0("Varimp for observation ", i, " saved."))
        }

        simes_pval.a <- simes.test(tau.a.stats[, 3])
        simes_pval.b <- simes.test(tau.b.stats[, 3])

        partial_simes_pval.a <- simes.partial(floor(no.obs.test * 0.05), tau.a.stats[, 3])
        partial_simes_pval.b <- simes.partial(floor(no.obs.train * 0.05), tau.b.stats[, 3])

        # check the correlation between two predictions from two datasets
        test_res.a <- correlation_test(tau.a.train.stats[, 1], tau.b.stats[, 1], methods = c('pearson', 'kendall', 'spearman'))
        test_res.b <- correlation_test(tau.a.stats[, 1], tau.b.train.stats[, 1], methods = c('pearson', 'kendall', 'spearman'))

        fisher.pval.a <- fisher.exact.test(tau.a.train.stats[, 3], tau.b.stats[, 3])
        fisher.pval.b <- fisher.exact.test(tau.a.stats[, 3], tau.b.train.stats[, 3])

        t.test.pval.a <- quantile.t.test(tau.a.train.stats[, 1], tau.b.stats[, 1])
        t.test.pval.b <- quantile.t.test(tau.a.stats[, 1], tau.b.train.stats[, 1]) 
        
        correlation_rslt <- rbind(c(simes_pval.a, partial_simes_pval.a, test_res.a, fisher.pval.a, t.test.pval.a), 
                                  c(simes_pval.b, partial_simes_pval.b, test_res.b, fisher.pval.b, t.test.pval.b))
        correlation_matrix = rbind(correlation_matrix, correlation_rslt)
    }


}