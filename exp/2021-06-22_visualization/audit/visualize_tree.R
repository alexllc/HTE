# script to visualize trees for treatment of interest

setwd("~/project/HTE/")

## source bin scripts
bin_ls = list.files("./bin")
for (bin in bin_ls){
    source(paste0("./bin/", bin))
}

library(treeshap)

# requires you to unload devtools (using ellipsis)?

project <- "KIRC"

if( !file.exists(paste0("./dat/wds_backup/", project, "_wds.csv")) ) {
    # prepare clinical and expression data for forest building
    clin <- fetch_clinical_data(cancer_type = project)

    # processing expression data just as the previous, except for scale, we want to see the trees
    exp <- fetch_exp_data(cancer_type = project, 
                        addBatch = T, 
                        numericBatch = T, 
                        scale = T, 
                        primaryTumorOnly = T, 
                        keepAllAliquot = F, 
                        formatPatient = T)
    exp$donorId <- rownames(exp)

    wds_scaled = inner_join(clin, exp , by = "donorId")
    wds_scaled = dplyr::select(wds_scaled, all_of(c("donorId","outcome", "age_at_initial_pathologic_diagnosis", "gender", "ajcc_pathologic_tumor_stage", "TSS", "portion", "plate", "center", tx_vector)))

    wds_scaled <- read.csv(paste0("./dat/wds_backup/", project, "_wds.csv"))

    cm_scaled= dplyr::select(wds_scaled, -c("donorId", "outcome"))

    write.csv(wds_scaled, paste0("./dat/wds_backup/", project, "_wds.csv"), row.names = F)
} else {
    wds_scaled <- read.csv(paste0("./dat/wds_backup/", project, "_wds_scaled.csv"))
}
# using the scaled versions
cm_scaled <- dplyr::select(wds_scaled, -c("donorId", "outcome"))

## DEA re-run 07-06-2020
DEGs = read.csv(paste0("./dat/tables/", project, "_DEGtable.csv"))
DEG_ls = as.character(DEGs$X)
tx_vector = DEG_ls[DEG_ls %in% colnames(exp)]
message(paste0("Treatments to assess: ", length(tx_vector)))

# Building the forest

# define parameters as before
obsNumber <- dim(cm_scaled)[1]
trainId <- sample(1: obsNumber, floor(obsNumber/2), replace = F)
registerDoParallel(10)

# CF construction
no.obs <- dim(cm_scaled)[1]
no.obs.train <- length(trainId)
no.obs.test <- no.obs - no.obs.train

Y <- wds_scaled$outcome

# Load treatment of interest based on our HTE expression results
sig_TE <- read.csv(paste0("/home/alex/project/HTE/exp/2021-02-03_proj_sum/res/2020-06-20_TCGA_pancancer_DEA_indicated_tx_dirct_.25_gp_as_tx_exp_HTE_", project, "_summary.csv"))
sig_TE <- filter(sig_TE, diff.pred.pval < 0.05 & mean.pred.pval < 0.05)
sig_TE_gene <- sig_TE$gene_names


# Build CF for each sig TE

best_tree_ls <- data.frame(gene = character(),
                           best_tree = double())
for (tx in sig_TE_gene) {
    cat(paste0(rep("=", 80), collapse = ""))
    cat(paste0("Processing: ", tx))
    if (sig_TE$logFC[which(sig_TE$gene_names == tx)] > 0) {
        treatment <- as.numeric(wds_scaled[[tx]] > quantile(wds_scaled[[tx]], 0.75))
    } else {
        treatment <- as.numeric(wds_scaled[[tx]] < quantile(wds_scaled[[tx]], 0.25))
    }

    X <- as.matrix(dplyr::select(cm_scaled, -all_of(tx)))

    CF <- cf(covarites = X, Y = Y, W = treatment, num_trees = 10000)   # run causal forests by default
    save(CF, file = paste0("/home/alex/project/HTE/exp/2021-06-22_visualization/dat/CF_", tx, ".RData"))
    # Verify sig TE
    tau <- predict(CF)
    tc_res <- test_calibration(CF)

    # Find best trees
    best_tree <- find_best_tree(CF, type = "causal")
    best_tree_ls <- append(best_tree_ls, c(list(tx), list(best_tree$best_tree)))

}

write.csv("")


# Calculate Shap values
unified <- randomForest.unify(rf_model = CF, data = X, W = treatment, Y = Y, is_grf = T)
