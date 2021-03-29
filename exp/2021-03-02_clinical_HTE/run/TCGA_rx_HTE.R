# Script for running HTE analysis using actual drug treatment as treatmetn variables

## Run from base directory
setwd("~/project/HTE/")

## source bin scripts
bin_ls = list.files("./bin")
for (bin in bin_ls){
    source(paste0("./bin/", bin))
}

## Set parameters for this run
cancerList <- c("BLCA", "COAD", "GBM", "HNSC", "KIRC", "LUAD","LUSC", "PRAD", "STAD", "THCA", "UCEC")

endpt <- "OS"
paused_at <- NULL
useALDEX2_DE <- TRUE
use_DE_covar_only <- TRUE
# Control strigency indicator
pureCtrlsOnly <- FALSE
save_wds <- TRUE

#****** Load various drug treatment information ******

for (cancer_type in cancerList) {
    message(rep("=", 50))
    message(paste0("Running cancer type: ", cancer_type))
    output_file <- paste0("./exp/2021-03-02_clinical_HTE/res/", cancer_type, "/")
    ## BCR biotab
    # Load drug taken from TCGA BCR Biotab files
    query <- GDCquery(project = paste0("TCGA-", cancer_type), 
                    data.category = "Clinical",
                    data.type = "Clinical Supplement", 
                    data.format = "BCR Biotab")
    GDCdownload(query)
    clinical.BCRtab.all <- GDCprepare(query)
    names(clinical.BCRtab.all)
    drug <- clinical.BCRtab.all[[paste0("clinical_drug_", tolower(cancer_type))]]


    ## Clinical indexed data
    clinical <- GDCquery_clinic(project = paste0("TCGA-", cancer_type), type = "clinical")

    ## surveying the dataset
    table(clinical$treatments_pharmaceutical_treatment_or_therapy)

        #   no not reported          yes 
        #  148          103          846

    nr = filter(clinical, treatments_pharmaceutical_treatment_or_therapy == "not reported")
    nr_drug = drug[drug$bcr_patient_barcode %in% nr$bcr_patient_barcode,]
    dim(nr)
        # [1] 103  74
    dim(nr_drug)
        # [1] 16 28
    length(unique(nr_drug$bcr_patient_barcode)) 
        # 10 patient has drug record despite being labelled as "not reported"


    no = filter(clinical, treatments_pharmaceutical_treatment_or_therapy == "no")
    no_drug = drug[drug$bcr_patient_barcode %in% no$bcr_patient_barcode,]
    dim(no)
    # [1] 148  74
    dim(nr_drug)
    # [1] 16 28
    length(unique(no_drug$bcr_patient_barcode)) 
        # 5 patient has drug record despite being labelled as "no"
        # 143 patients are truly not taking drugs

    yes = filter(clinical, treatments_pharmaceutical_treatment_or_therapy == "yes")
    yes_drug = drug[drug$bcr_patient_barcode %in% yes$bcr_patient_barcode,]
    dim(yes)
    # [1] 846  74
    dim(yes_drug)
    # [1] 2378   28
    length(unique(yes_drug$bcr_patient_barcode))
        # 765 
        # 99 patients don't have drug names recorded despite being labelled as "yes"

    # No: 143
    # Not reported: 93 + 99 = 192
    # Yes: 10 + 5 + 765 = 780

    # For each drug, we can choose the truely "no" patients as control, these are the universally non-treated patiets. If a patient with drug recrods is not treated with the drug in question but treated with other drugs in question, they will be the pseudo control group. We could try HTE with pseudo control or true control, depending on whichever works.

    # Fetch clinical indexed data from TCGAbiolinks
    clin_indexed <- fetch_indexed_clinical(cancerType = cancer_type)
    clin_indexed <- convert_col_to_numeric(clin_indexed, id = "bcr_patient_barcode")

    # Now the next step is to set up the drug binary tables, one hot encoding based on the above criteria

    drug_taken <- dplyr::select(drug, c("bcr_patient_barcode", "pharmaceutical_therapy_drug_name", "pharmaceutical_therapy_type"))
    drug_taken <- drug_taken[-c(1:2),]

    # Unify drug names with pre-downloaded conversion table
    drug_taken$corr_drug_names <- correct_drug_names(drug_taken$pharmaceutical_therapy_drug_name)
    drug_taken$corr_drug_names <- gsub(" ", "", drug_taken$corr_drug_names)
    drug_taken <- dplyr::select(drug_taken, -"pharmaceutical_therapy_drug_name") # get the old names out of the way

    ## Universal controls
    controls <- no$bcr_patient_barcode
    controls <- unique(controls[!(controls %in% drug_taken$bcr_patient_barcode)])

    ## True not reported cases (exclude)
    not_reported <- nr[!(nr$bcr_patient_barcode %in% drug_taken$bcr_patient_barcode),]
    not_reported <- not_reported$bcr_patient_barcode

    ## Load other covariates for covarriate matrix X input
    iqlr_count <- fread(paste0("./dat/iqlr_normalized_expr/TCGA-", cancer_type, "_iqlr_expected_count.csv.gz"))

    # Differentially expressed genes included in the covar only
    de <- read.csv(paste0("./dat/tables/", cancer_type, "_DEGtable.csv"))
    iqlr_count <- iqlr_count[iqlr_count$V1 %in% de$X,]
    colnames(iqlr_count)[colnames(iqlr_count) == "V1"] <- "gene"
    colnames(iqlr_count) <- gsub("\\.", "-", colnames(iqlr_count))
    selected_aliq <- filter_replicate_samples(colnames(iqlr_count)[2:dim(iqlr_count)[2]])
    iqlr_count <- dplyr::select(iqlr_count, all_of(c(c("gene", selected_aliq))))
    iqlr_count <- t(iqlr_count)
    colnames(iqlr_count) <- iqlr_count[1,]
    iqlr_count <- iqlr_count[-1,]
    mode(iqlr_count) = "numeric"
    iqlr_count <- as.data.frame(iqlr_count)
    iqlr_count$donorId <- format_tcga_patient(rownames(iqlr_count))
    head(iqlr_count)[,1:10]

    # head(sig)[, c(1:3, 1219:1226)]

    # co-prescription one hot encoding
    hot_drug <- dplyr::select(drug_taken, -c("pharmaceutical_therapy_type"))

    # re-enter multi-drug treated patients into multiple rows
    multi_tx <- hot_drug[grep("\\+", hot_drug$corr_drug_names),]

    for (i in 1:nrow(multi_tx)) {
        record <- unlist(strsplit(multi_tx$corr_drug_names[i], "\\+"))
        for (j in 1:length(record)) {
            hot_drug <- rbindlist(list(hot_drug, as.list(c(multi_tx$bcr_patient_barcode[i], record[j]))))
        }
    }
    hot_drug <- data.table(hot_drug[!grep("\\+", hot_drug$corr_drug_names)]) # removing rows by index not supported in data.tables
    hot_drug <- dcast(setDT(hot_drug), bcr_patient_barcode ~ corr_drug_names, value.var = "corr_drug_names", function(x) as.numeric(length(x) > 0))

    # add true controls to same tabel
    control_row <- unique(no_drug$bcr_patient_barcode)

    for (ctrl in control_row) {
        hot_drug <- rbindlist(list(hot_drug, c(ctrl, as.list(rep(0, dim(hot_drug)[2] - 1)))))
    }

    # Bind potentially co-administered drug to clinical table
    clin_indexed <- inner_join(clin_indexed, hot_drug, by = "bcr_patient_barcode")


    # result saving headers
    tc_res <- NULL
    ape_res <- NULL
    blp_res <- NULL
    ate_res <- NULL

    if (!is.null(paused_at)) {
        drug_selection <- colnames(hot_drug)[colnames(hot_drug) != "bcr_patient_barcode"][1:(paused_at - 1)] # to resume drug from previous run
        # drug_selection <- colnames(hot_drug)[colnames(hot_drug) != "bcr_patient_barcode"][(paused_at - 1):(length(colnames(hot_drug)) - 1)]
    } else {
        drug_selection <- colnames(hot_drug)[colnames(hot_drug) != "bcr_patient_barcode"]
    }

    for (rx in  drug_selection) {

        message(paste0(rep("=", 50)))

        message(paste0("Running drug: ", rx))
        drug_takers <- unique(drug_taken$bcr_patient_barcode[drug_taken$corr_drug_names == rx])
        # a control for the current rx, but these patients take other drugs
        pseudo_ctr <- unique(drug_taken$bcr_patient_barcode[!(drug_taken$bcr_patient_barcode %in% not_reported)])

        wds <- left_join(clin_indexed, iqlr_count, by = c("bcr_patient_barcode" = "donorId"))
        if (save_wds) {
            saveRDS(wds, file = paste0(output_file, cancer_type, "_wds.rds.gz"))
            message("Whole dataset prepared and saved.")
        }
        treatment_W <- wds[[rx]]
        print(table(treatment_W))

        X <- dplyr::select(wds, -c(rx, "bcr_patient_barcode", "OS_time", "vital_status"))

        # we're not sure if each drug had enough observations, so sometimes this might fail upon forest builing
        try(forest <- causal_forest(X = X, 
                                Y = wds$OS_time, 
                                W = treatment_W, 
                                num.trees = 20000, 
                                seed = 2020,
                                tune.parameters = c("sample.fraction", "mtry", "min.node.size", "honesty.fraction", "honesty.prune.leaves", "alpha", "imbalance.penalty")
                            )
            )
        if(!exists("forest")) {
            message("Failed to build forest.")
            next
        }

        # cf_scores <- grf::get_scores(cf) # does not work
        predictions <- predict(forest, estimate.variance = TRUE)

        # try(cf_split_freq <- split_frequencies(forest, max.depth = 4))
        try(cf_tc <- test_calibration(forest))
        try(cf_ape <- average_partial_effect(forest))
        try(cf_blp <- best_linear_projection(forest))
        try(cf_ate <- average_treatment_effect(forest))
        
        if( !all(c(exists("cf_tc"),exists("cf_ape"),exists("cf_blp"),exists("cf_ate"))) | any(c( any(is.na(cf_tc)), any(is.na(cf_ape)), any(is.na(cf_blp)), any(is.na(cf_ate)) ) ) ) {
            message("Failed to analyse forest.")
            next
        }

        message("Excess error summary statistics:")
        print(summary(predictions$excess.error))
        print(cf_tc)
        print(cf_ape)
        print(cf_blp)
        print(cf_ate)

        # save res for rbindlist later
        tc_res <- rbind(tc_res, c(rx, cf_tc[1,], cf_tc[2,]))
        ape_res <- rbind(ape_res, c(rx,cf_ape))
        ate_res <- rbind(ate_res, c(rx,cf_ate))
        blp_res <- rbind(blp_res, c(rx,cf_blp)) # estimation of beta0 only

        # estimate feature importance
        varimp <- variable_importance(forest)
        varimp <- as.data.frame(cbind(colnames(X),varimp))
        colnames(varimp) <- c("feature", "importance")
        varimp <- varimp[order(varimp$importance, decreasing = TRUE),]

        # translate gene names in varimp into readable form
        important_gene <- varimp[grep("ENSG*", varimp$feature),1] # only keep genes
        # Convert ensemblIDs to ENTERZ gene ID for the query
        suppressMessages(ensembl2ID <- AnnotationDbi::select(org.Hs.eg.db, keys = important_gene,columns=c("SYMBOL"), keytype="ENSEMBL"))

        varimp <- left_join(varimp, ensembl2ID, by = c("feature" = "ENSEMBL"))

        write.csv(varimp, file = paste0(output_file, rx, "_varimp.csv"), row.names = FALSE)

        # Calculate BLP while considering top 1% of important features
        important_feature = varimp$feature[as.numeric(varimp$importance) > quantile(as.numeric(varimp$importance), 0.99)]
        cf_blp_A <- try(best_linear_projection(forest, A = X[, important_feature]))
        if (class(cf_blp_A) == "try-error") next
        write.csv(cf_blp_A, file = paste0(output_file, rx, "_BLP_Ai.csv"))

    }
    # Save pan-drug results
    write.xlsx(tc_res, file = paste0(output_file, cancer_type, "_grf_res.xlsx"), sheetName="test_calib", row.names=FALSE)
    write.xlsx(ape_res, file=paste0(output_file, cancer_type, "_grf_res.xlsx"), sheetName="avg_partial_eff", append=TRUE, row.names=FALSE)
    write.xlsx(blp_res, file=paste0(output_file, cancer_type, "_grf_res.xlsx"), sheetName="best_linear_proj_beta_0", append=TRUE, row.names=FALSE)
    write.xlsx(ate_res, file=paste0(output_file, cancer_type, "_grf_res.xlsx"), sheetName="avg_tx_eff", append=TRUE, row.names=FALSE)

    message("Drug HTE analysis completed. Beginning SHC and permutation tests.")

    # Run Kai's HTE tests
    # no need to try here, if ti fails it will fail up top
    obsNumber <- dim(wds)[1]
    trainId <- sample(1: obsNumber, floor(obsNumber/2), replace = FALSE)
    registerDoParallel(10)
    colnames(wds)[colnames(wds) == "bcr_patient_barcode"] <- "donorId"
    colnames(wds)[colnames(wds) == "OS_time"] <- "outcome"
    result <- run.hte(covar_mat = dplyr::select(wds, -c("donorId", "outcome", "vital_status")), 
                    tx_vector = colnames(hot_drug)[colnames(hot_drug) != "bcr_patient_barcode"], 
                    whole_dataset = wds, 
                    project = cancer_type, 
                    W_matrix = wds[,colnames(hot_drug)[colnames(hot_drug) != "bcr_patient_barcode"]], 
                    trainId = trainId, 
                    seed = 111, 
                    is_binary = T, 
                    is_save = T, 
                    save_split = T, 
                    is_tuned = F, 
                    thres = 0.75, 
                    n_core = 70, 
                    output_directory = paste0(output_file, "perm_shc/"), 
                    perm_all = TRUE,
                    random_rep_seed = TRUE
                    ) # pre-filtered thres, should not have an effect here
    write.csv(result[[1]], paste0(output_file, cancer_type, '_drug_correlation_test_result.csv'), quote = F, row.names = F)
    write.csv(result[[2]], paste0(output_file, cancer_type, '_drug_calibration_result.csv'), quote = F, row.names = F)
    write.csv(result[[3]], paste0(output_file, cancer_type, '_drug_median_t_test_result.csv'), quote = F, row.names = F)
    write.csv(result[[4]], paste0(output_file, cancer_type, '_drug_permutate_testing_result.csv'), quote = F, row.names = F)
    write.csv(result[[5]], paste0(output_file, cancer_type, '_drug_obs_tau_risk_var.csv'), quote = F, row.names = F)
    write.csv(result[[6]], paste0(output_file, cancer_type, '_drug_avg_tx_effect.csv'), quote = F, row.names = F)
    write.csv(result[[7]], paste0(output_file, cancer_type, '_drug_avg_partial_tx_effect.csv'), quote = F, row.names = F)
    write.csv(result[[8]], paste0(output_file, cancer_type, '_drug_best_linear_pred_intercept.csv'), quote = F, row.names = F)

}


# The TCGA-CDR summarized clinical data collected before the publication date. However, new data is being uploaded to the GDC data portal constantly. In this example, patient "TCGA-OL-A66H" was added to the database at 2019-08-08T16:30:14.878156-05:00, who would not be included in the TCGA-CDR of course. Therefore, for the missing patients, you would have to 'scavenge' the other datasets for their survival length.

    # Quoted from https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/clinical.html

    # In GDC database the clinical data can be retrieved from different sources:

    #     indexed clinical: a refined clinical data that is created using the XML files.
    #     XML files: original source of the data
    #     BCR Biotab: tsv files parsed from XML files

    # There are two main differences between the indexed clinical and XML files:

    #     XML has more information: radiation, drugs information, follow-ups, biospecimen, etc. So the indexed one is only a subset of the XML files
    #     The indexed data contains the updated data with the follow up information. For example: if the patient is alive in the first time clinical data was collect and the in the next follow-up he is dead, the indexed data will show dead. The XML will have two fields, one for the first time saying he is alive (in the clinical part) and the follow-up saying he is dead.


    # 1. Clinical variables
    # Suppose we take the clinical indexed file from faith (we will verify potential conflicts later)
    
    # 2. Other drugs taken

    # 3. Differentially expressed genes normalized by IQLR


