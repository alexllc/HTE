library(cgdsr)
# setwd("../../HTE")
project = "BRCA"
source("./Overlap_of_varImp_ver2.R")

aggr_res <- function(res_mat, i, est_col_list = NULL) {
    removed_na <- na.omit(res_mat[, i])
    # simes_pval <- ifelse(length(removed_na) > 0, simes.test(removed_na), NA)
    if(!(i %in% est_col_list)){
        pval <- ifelse(length(removed_na) > 0, simes.partial(2, removed_na), NA)
    } else {
        pval <- ifelse(length(removed_na) > 0, extract_binom_pval(removed_na), NA) 
    }
}

binary_tx <- function(treatment, tx_dirct, tx_gene, thres) {
    if(tx_dirct[tx_gene] > 0) {
        treatment = as.numeric(treatment > quantile(treatment, thres))
    } else {
        treatment  = as.numeric(treatment < quantile(treatment, 1- thres))}
}

gene2ensembl <- function(genels, from  = "ENSEMBL", to = "SYMBOL") {
    # use keytypes(org.Hs.eg.db) to list available conversion types
    txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
    txnames =AnnotationDbi::select(Homo.sapiens, keys = genels, columns = to, keytype = from, multiVals = "CharacterList")
    genels = as.data.frame(genels)
    genels = left_join(genels, txnames, by = c("genels" = from)) %>% group_by(genels) %>% dplyr::slice(1) %>% as.data.frame() # multiple matches are reduced by picking the first match
    return(as.character(genels[,to]))
}

format_tcga_patient <- function(pat_ls) {
    tmp = strsplit(pat_ls, "-")
    tmp = unlist(lapply(tmp, function(x) paste(x[[1]], x[[2]], x[[3]], sep = "-")))
    return(unlist(tmp))
}


mode_list = c("mutation", "microarray")
mode = mode_list[2]


if (mode == "mutation") {

    source("./HTE_mutation.R")
    output_file = "./result/BRCA/"

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
    whole_dataset = left_join(tcga_clin, dplyr::select(wtcga, all_of(c("donorId", sel_genes))), by = "donorId")
    whole_dataset = dplyr::select(whole_dataset, -c("gender", "ajcc_pathologic_tumor_stage"))
    whole_dataset = whole_dataset[complete.cases(whole_dataset),]
    covar_mat= dplyr::select(whole_dataset, -c("donorId", "outcome"))


    ## Running HTE for both datasets
    covar_type = "mutation"
} else if (mode == "microarray") {
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

    # For conversion between ensembl and gene symbol
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(Homo.sapiens)

    setwd("/home/alex/project/HTE/wd/HTE/")
    source("./grf_parameters.R")
    source("./HTE_main_functions.R")
    source("./HTE_validation_functions.R")
    source("./survival_imputation.R")
    setwd("/home/alex/project/HTE/wd/validation/METABRIC/METAB-TCGA-expression")

    output_file = "/home/alex/project/HTE/wd/expression_HTE/result/METABRIC_SHC_TCGA/"
    metab_dirt = "/home/alex/project/HTE/wd/mut_HTE/METABRIC/data_jun11/"
    tcga_dirt = "/home/alex/project/HTE/wd/validation/cbioportal_tcga_brca_pancan/"
    
    clin_col = c("donorId", "age", "outcome")

    ## (1) a. Get and process METABRIC clinical data. Same processing step as the METABRIC_microarray.R file, please refer to original for documentation
    samp = fread(paste0(metab_dirt, 'data_clinical_sample.txt'), skip=4)
    samp = samp[,c(1,4:7,9)]
    pat = fread(paste0(metab_dirt, 'data_clinical_patient.txt'), skip=4)
    pat = pat[,c(1:3,5,7:10,12,13,14,19)]
    pat = pat[which( !(pat$OS_MONTHS <= 0 | is.na(pat$OS_MONTHS) | pat$OS_STATUS == "")),]
    pat$OS_STATUS = as.integer(pat$OS_STATUS=="1:DECEASED") # METABRIC updated the OS status term
    pat$OS_MONTHS = pat$OS_MONTHS * 30.4167
    attach(pat)
    imp.dat = NNMIS(NPI, xa = AGE_AT_DIAGNOSIS, xb = AGE_AT_DIAGNOSIS, time = OS_MONTHS, event = OS_STATUS, imputeCT = T, Seed = 2020, mc.cores = 30)
    detach(pat)
    imp_surv = imp.dat$dat.T.NNMI %>% mutate(mean = rowMeans(.))
    pat$outcome = imp_surv$mean
    pat$OS_MONTHS = NULL
    pat$OS_STATUS = NULL
    metab_clin = merge(pat, samp, by = "PATIENT_ID")
    metab_clin = as.data.frame(metab_clin[complete.cases(metab_clin),])
    for (c in colnames(metab_clin)) {

        if (!is.numeric(metab_clin[,c]) && c != "PATIENT_ID") {
            which.one <- which( levels(metab_clin[,c]) == "")
            levels(metab_clin[,c])[which.one] <- NA
            metab_clin[,c] = sapply(sapply(metab_clin[,c], as.factor), as.numeric) 
            print(paste0(c, " is converted to numeric.")) 
        }
    }
    metab_clin = dplyr::select(metab_clin, PATIENT_ID, AGE_AT_DIAGNOSIS, outcome)
    colnames(metab_clin) = clin_col

    ## (1) b. Import, transpose and save z-score data for METABRIC patients. Detail about normalization process is found in: https://docs.cbioportal.org/5.1-data-loading/data-loading/file-formats/z-score-normalization-script. We're using z-scores with ref. to diploid samples.
    if(!file.exists("./metab_transposed_data_expression_median_zscore.txt")) {
        metab_z = fread(paste0(metab_dirt,"data_mRNA_median_Zscores.txt"))
        message("Transposing expression z-score matrix, this could take a *LONG* while.")
        metab_z = t(metab_z)
        colnames(metab_z) = metab_z[1,]
        metab_z = metab_z[-c(1:2),]
        class(metab_z) = "numeric"
        metab_z = as.data.frame(metab_z)
        metab_z$donorId = rownames(metab_z)
        write.table(metab_z, "metab_transposed_data_expression_median_zscore.txt", sep = "\t", row.names = F) # as of R4.0 R no longer automatically converts strings as factors
    } else {
        metab_z = as.data.frame(fread("metab_transposed_data_expression_median_zscore.txt"))
    }
    ## METABRIC data sets does not have fully NA columns
    metab_not_na = sapply(metab_z, function(x)all(!is.na(x)))
    metab_z = metab_z[,metab_not_na] # only 1013 genes that are complete

    ## (1) c. Import metab_txdirct infomraiton about METABRIC expression for treatment direciton
    metab_txdirct = fread(paste0(metab_dirt, "data_CNA.txt"))
    metab_txdirct = as.data.frame(t(as.matrix(metab_txdirct)))
    colnames(metab_txdirct) = metab_txdirct[1,]
    metab_txdirct = metab_txdirct[-c(1,2),]
    metab_txdirct = metab_txdirct[complete.cases(metab_txdirct),]
    metab_txdirct = sapply(metab_txdirct, as.numeric)
    metab_txdirct = apply(metab_txdirct, 2, sum)
    names(metab_txdirct) = gene2ensembl(names(metab_txdirct), from = "SYMBOL", to = "ENSEMBL")

    ## (2) a. Get and process TCGA clinical data
    tcga_clin = read.csv("/home/alex/project/HTE/wd/expression_HTE/TCGA_CDR_clean.csv")
    tcga_clin = filter(tcga_clin, type == "BRCA")
    labels = c("[Discrepancy]","[Not Applicable]","[Not Available]","[Unknown]")
    tcga_clin$ajcc_pathologic_tumor_stage[which(tcga_clin$ajcc_pathologic_tumor_stage %in% labels)] = NA

    for (c in colnames(tcga_clin)) {
        if (!is.numeric(tcga_clin[,c]) && c != "donorId") {
            which.one <- which( levels(tcga_clin[,c]) == "")
            levels(tcga_clin[,c])[which.one] <- NA
            tcga_clin[,c] = sapply(sapply(tcga_clin[,c], as.factor), as.numeric) 
            print(paste0(c, " is altered")) 
        }
    }
    attach(tcga_clin)
    tcga_imp = NNMIS(ajcc_pathologic_tumor_stage, xa = age_at_initial_pathologic_diagnosis, xb = age_at_initial_pathologic_diagnosis, time = OS.time, event = OS, imputeCT = T, Seed = 2020, mc.cores = 30)
    detach(tcga_clin)
    tcga_imp_surv = tcga_imp$dat.T.NNMI %>% mutate(mean = rowMeans(.))
    tcga_clin$outcome = tcga_imp_surv$mean
    tcga_clin = dplyr::select(tcga_clin, donorId, age_at_initial_pathologic_diagnosis, outcome)
    colnames(tcga_clin) = clin_col

    ## (2) b. Import z-score for TCGA patients downlaoded from cbioportal pancancer atlas.
    if(!file.exists("./tcga_transposed_data_expression_median_zscore.txt")) {
        tcga_z = fread(paste0(tcga_dirt,"data_RNA_Seq_v2_mRNA_median_Zscores.txt"))
        message("Transposing expression z-score matrix, this could take a *LONG* while.")
        tcga_z = t(tcga_z)
        colnames(tcga_z) = tcga_z[1,]
        tcga_z = tcga_z[-c(1:2),]
        class(tcga_z) = "numeric"
        tcga_z = as.data.frame(tcga_z)
        tcga_z$donorId = rownames(tcga_z)
        write.table(tcga_z, "tcga_transposed_data_expression_median_zscore.txt", sep = "\t", row.names = F) # as of R4.0 R no longer automatically converts strings as factors
    } else {
        tcga_z = as.data.frame(fread("tcga_transposed_data_expression_median_zscore.txt")) # donorId will be right at the back
    }
    tcga_not_na = sapply(tcga_z, function(x)all(!is.na(x))) # many genes are fully NA
    tcga_z = tcga_z[,tcga_not_na]
    tcga_z$donorId = format_tcga_patient(tcga_z$donorId)
    ## All are complete cases
    # tcga_z = tcga_z[complete.cases(tcga_z),]

    ## (2) c. Import DEA results for treatment direction
    T_DEG = read.csv(paste0("/home/alex/project/HTE/wd/expression_HTE/tables/", project, "_DEGtable.csv"))
    T_txdirct = T_DEG$logFC
    names(T_txdirct) = gene2ensembl(T_DEG$X, from = "ENSEMBL", to = "SYMBOL")

    # Select only the overlapping genes
    overlap_genes = intersect(colnames(metab_z), colnames(tcga_z))
    overlap_genes = overlap_genes[-which(overlap_genes == "donorId")]
    
    metab_whole = dplyr::select(metab_z, all_of(c("donorId", overlap_genes))) %>% left_join(metab_clin, metab_whole, by = "donorId")
    metab_whole = metab_whole[complete.cases(metab_whole),]
    tcga_whole = dplyr::select(tcga_z, all_of(c("donorId", overlap_genes))) %>% left_join(tcga_clin, by = "donorId")
    tcga_whole = tcga_whole[complete.cases(tcga_whole)]
    metab_whole = dplyr::select(metab_whole, all_of(c("donorId", "outcome", overlap_genes)))
    tcga_whole = dplyr::select(metab_whole, all_of(c("donorId", "outcome", overlap_genes)))
    metab_covar = dplyr::select(metab_whole, -c("donorId", "outcome"))
    tcga_covar= dplyr::select(tcga_whole, -c("donorId", "outcome"))

    ## Select significant treatment genes 
    overlap_DEG = intersect(names(T_txdirct), names(metab_txdirct)) # 166 for BRCA
    overlap_DEG = overlap_DEG[overlap_DEG %in% overlap_genes]
}

## SHC code for the prepared data
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
overlap_test_names = c("gene","fisher_top1_pval","fisher_top3_pval","fisher_top5_pval","fisher_top10_pval","fisher_top20_pval","fisher_top30_pval","above0_pval","perm_fisher_pval","perm_min_pval","perm_simes_pval","perm_softomni_pval")

for (tx in overlap_DEG){
    
    print(paste0(c('#', rep('-', 40), ' running ', which(overlap_DEG == tx), ' of ', length(overlap_DEG), rep('-', 40)), collapse = ''))
    # tx vector and covar mat for TCGA
    cat('Treatment name:', tx, fill = T)

    T_treatment <- as.data.frame(tcga_covar)[, tx]
    T_covariates <- as.matrix(dplyr::select(tcga_covar, -tx))
    ## CHANGE WHEN RUNNING EXPRESSION/MUTATION
    if (mode == "mutation") {
        T_treatment <- as.numeric(T_treatment != 0) # only for mutation
    } else {
       T_treatment <- binary_tx(treatment = T_treatment, tx_dirct = T_txdirct, tx_gene = tx, thres = thres)
    }
    T_Y <- tcga_whole$outcome

    # tx vector and covar mat for METAB
    M_treatment <- as.data.frame(metab_covar)[, tx]
    M_covariates <- as.matrix(dplyr::select(metab_covar, -tx))
    
    if (mode == "mutation") {
        M_treatment <- as.numeric(M_treatment != 0) # only for mutation
    } else {
       M_treatment <- binary_tx(treatment = M_treatment, tx_dirct = metab_txdirct, tx_gene = tx, thres = thres)
    }
    
    M_Y <- metab_whole$outcome

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

    # Skipping repeats, seemes unnecessary
    # no_repeats = 10
    # for (i in seq(no_repeats)) {
    i = 0
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

    if( all(test_res_tcga[c(2,4,6)] < 0.05)) {
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
        write.csv(overlap_matrix, file = paste0(file_prefix, '_varimp_overlap.csv'), row.names = F, quote = F)
        aggregated_overlap_rslt <- sapply(seq(dim(overlap_matrix)[2]), aggr_res, est_col_list = 7, res_mat = overlap_matrix)
        message("===========Varimp overlap results: ==============")
        for(k in 2:length(aggregated_overlap_rslt)) {
            cat(paste0(overlap_test_names[k], ": "), aggregated_overlap_rslt[k], fill = T)
        }
    } else {
        message("TCGA SHC test correlation not significant, skipping varimp overlap test.")
    }

    if(is_save){
        colnames(correlation_matrix) <- col_names 
        write.csv(correlation_matrix, file = paste0(file_prefix, '_split_half.csv'), row.names = F, quote = F)
    }

    aggregated_corr_rslt <- sapply(seq(dim(correlation_matrix)[2]), aggr_res, est_col_list = c(3, 5, 7), res_mat = correlation_matrix)


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


    # Alex @ Jun28, 2020 do.call then append to preset data.frame does not work anymore, switching to rbindlist
    current_ret <- do.call('c', list(list(tx), as.list(aggregated_corr_rslt[1: 8])))
    current_ret <- rbindlist(list(current_ret))
    correlation_test_ret <- rbind(correlation_test_ret, current_ret)

    overlap_tmp <- do.call('c', list(list(tx), as.list(aggregated_overlap_rslt)))
    overlap_tmp <- rbindlist(list(overlap_tmp))
    overlap_test_res <- rbind(overlap_test_res, overlap_tmp)
                             
}

colnames(correlation_test_ret) <- corr_test_names
colnames(overlap_test_res) <- overlap_test_names
write.csv(correlation_test_ret, paste0(output_file, project, '_validation_correlation_test_result.csv'), quote = F, row.names = F)
write.csv(overlap_test_res, paste0(output_file, project, '_validation_varimp_overlap_result.csv'), quote = F, row.names = F)