# Script testing whether covariates are stabaliy selected over different seeds

## Run from base directory
setwd("~/project/HTE/")

## source bin scripts
bin_ls = list.files("./bin")
for (bin in bin_ls){
    source(paste0("./bin/", bin))
}

## Set parameters for this run "BLCA", "COAD",  

endpt <- "OS"
paused_at <- NULL
use_DE_covar_only <- TRUE
# Control strigency indicator
pureCtrlsOnly <- FALSE
save_wds <- TRUE

cancer_type <- "BRCA"

hot_drug <- fetch_drug(cancer_type)

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

for (rx in  drug_selection) {

    message(paste0(rep("=", 50)))

    message(paste0("Running drug: ", rx))
    drug_takers <- unique(drug_taken$bcr_patient_barcode[drug_taken$corr_drug_names == rx])
    # a control for the current rx, but these patients take other drugs
    pseudo_ctr <- unique(drug_taken$bcr_patient_barcode[!(drug_taken$bcr_patient_barcode %in% not_reported)])

    wds <- left_join(clin_indexed, iqlr_count, by = c("bcr_patient_barcode" = "donorId"))
    treatment_W <- wds[[rx]]
    print(table(treatment_W))
    if (all(!as.logical(treatment_W))) {
        message("No drug takers with both drug and clinical information")
        next
    }
    if (save_wds) {
        saveRDS(wds, file = paste0(output_file, cancer_type, "_wds.rds.gz"))
        message("Whole dataset prepared and saved.")
    }

    X <- dplyr::select(wds, -c(rx, "bcr_patient_barcode", "OS_time", "vital_status"))

    # we're not sure if each drug had enough observations, so sometimes this might fail upon forest builing
    varimp_99_record <- 
    for (k in 1:num_rep) {
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

        # estimate feature importance
        varimp <- variable_importance(forest)
        varimp <- as.data.frame(cbind(colnames(X),varimp))
        colnames(varimp) <- c("feature", "importance")
        varimp <- varimp[order(varimp$importance, decreasing = TRUE),]
        important_feature = varimp$feature[as.numeric(varimp$importance) > quantile(as.numeric(varimp$importance), 0.99)]
    }