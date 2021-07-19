# script to visualize trees for treatment of interest

setwd("~/project/HTE/")

## source bin scripts
bin_ls = list.files("./bin")
for (bin in bin_ls){
    source(paste0("./bin/", bin))
}

library(treeshap)
library(compositions)

# requires you to unload devtools (using ellipsis)?

project <- "KIRC"

if( !file.exists(paste0("./dat/wds_backup/", project, "_wds.csv")) ) {
    # prepare clinical and expression data for forest building
    clin <- fetch_clinical_data(cancer_type = project)

    ## DEA re-run 07-06-2020
    DEGs = read.csv(paste0("./dat/tables/", project, "_DEGtable.csv"))
    DEG_ls = as.character(DEGs$X)

    # using IQLR transformed expression data instead of FPKM-UQ for post-hoc analysis
    iqlr_exp <- fread(paste0("dat/iqlr_normalized_expr/", "TCGA-", project, "_iqlr_expected_count.csv.gz"))
    patient_id <- colnames(iqlr_exp)
    gene_id <- iqlr_exp$V1

    # preselection makes the transpose much faster
    sel_DEG_ls <- DEG_ls[DEG_ls %in% gene_id] 
    setkey(iqlr_exp, V1)
    iqlr_exp <- iqlr_exp[sel_DEG_ls]
    gene_id <- iqlr_exp$V1 # update gene id list
    
    iqlr_exp <- iqlr_exp[,-1] # remove the gene names
    iqlr_exp <- t(iqlr_exp)
    colnames(iqlr_exp) <- gene_id

    # format and select patient ID
    rownames(iqlr_exp) <- gsub("\\.", "-", rownames(iqlr_exp))
    iqlr_exp <- iqlr_exp[filter_replicate_samples(rownames(iqlr_exp)),]

    # add batch as covariates
    iqlr_exp <- cbind(separate(as.data.frame(rownames(iqlr_exp)), 
                            "rownames(iqlr_exp)", 
                            c(NA, "TSS", NA, NA, "portion", "plate", "center"),  # skip var with NAs
                            sep = "-"), 
                    iqlr_exp)

    # convert batches to numeric values
    batches = c("TSS", "patient", "portion", "plate", "center")
    for (name in batches) {
        c = grep(name, colnames(iqlr_exp))
        which.one <- which( levels(iqlr_exp[,c]) == "")
        levels(iqlr_exp[,c])[which.one] <- NA
        print(paste0(colnames(iqlr_exp)[c], " is converted to numeric"))
        iqlr_exp[,c] = sapply(sapply(iqlr_exp[,c], as.factor), as.numeric) 
    } 

    donorId <- format_tcga_patient(rownames(iqlr_exp))
    iqlr_exp <- cbind(donorId, iqlr_exp)
    
    # combine to wds
    wds <- left_join(clin, iqlr_exp, by = "donorId")
    covar <- dplyr::select(wds, -c("donorId", "outcome"))

    # missing data handling


    write.csv(wds, file = gzfile(paste0("./dat/wds_backup/", project, "_iqlr_wds.csv.gz")), row.names = F)
} else {
    wds <- read.csv(paste0("./dat/wds_backup/", project, "_iqlr_wds.csv.gz"))
}
tx_vector <- DEG_ls[DEG_ls %in% colnames(wds)]
message(paste0("Treatments to assess: ", length(tx_vector)))


# Building the forest

# define parameters as before
obsNumber <- dim(covar)[1]
trainId <- sample(1: obsNumber, floor(obsNumber/2), replace = F)
registerDoParallel(10)

# CF construction
no.obs <- dim(covar)[1]
no.obs.train <- length(trainId)
no.obs.test <- no.obs - no.obs.train

Y <- wds$outcome

# Load treatment of interest based on our HTE expression results
sig_TE <- read.csv(paste0("/home/alex/project/HTE/exp/2021-02-03_proj_sum/res/2020-06-20_TCGA_pancancer_DEA_indicated_tx_dirct_.25_gp_as_tx_exp_HTE_", project, "_summary.csv"))
sig_TE <- filter(sig_TE, diff.pred.pval < 0.05 & mean.pred.pval < 0.05)
sig_TE_gene <- sig_TE$gene_names

# check if missing at random
missing_sum <- missingSummary(wds)

# Build CF for each sig TE

best_tree_ls <- data.frame(gene = character(),
                           best_tree = double())
for (tx in sig_TE_gene) {
    cat(paste0(rep("=", 80), collapse = ""))
    cat(paste0("Processing: ", tx))

    # handling missing values


    if (sig_TE$logFC[which(sig_TE$gene_names == tx)] > 0) {
        treatment <- as.numeric(wds[[tx]] > quantile(wds[[tx]], 0.75))
    } else {
        treatment <- as.numeric(wds[[tx]] < quantile(wds[[tx]], 0.25))
    }

    X <- as.matrix(dplyr::select(covar, -all_of(tx)))

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
