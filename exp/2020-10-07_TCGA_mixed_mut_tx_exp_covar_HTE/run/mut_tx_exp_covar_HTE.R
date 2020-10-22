#' Script to run an HTE with mutation as treatment variable and expression as running covariates
# Run from base directory
setwd("~/project/HTE/")

library(metap)

# source bin scripts
bin_ls = list.files("./bin")
for (bin in bin_ls){
    source(paste0("./bin/", bin))
}

# Set param for this run
cancer_type = "BRCA"

# Retreive data
exp = fetch_exp_data(cancer_type, scale = FALSE)
exp = mk_id_rownames(exp)

head(exp)[,1:10]

endpt_ls = c("OS", "DSS", "DFI", "PFI")

if(!file.exists(paste0("./dat/mutation_frequencies/MAF_wide_", cancer_type, ".csv.gz"))) {
    mut = as.data.frame(fetch_mut_data(cancer_type)) # will take a long time, you'd want to save it
    write.csv(mut, file = gzfile(paste0("./dat/mutation_frequencies/MAF_wide_", cancer_type, ".csv.gz"), row.names = FALSE))
} else {
    mut = as.data.frame(fread(paste0("./dat/mutation_frequencies/MAF_wide_", cancer_type, ".csv.gz")))
}
mut = mk_id_rownames(mut)

head(mut)[,1:10]

set.seed(2020)

for (endpt in endpt_ls) {
    message(paste0("Using ", endpt, " as endpoint measurement.")
    clinical = fetch_clinical_data(cancer_type, outParam = endpt)
    clinical = mk_id_rownames(clinical)
    head(clinical)[,1:10]
    assign(paste0("clinical_", endpt), clinical) # saving clinical df for each endpoint types rather than a vector, might change to numeric vectors later
    
    # Check the list of common patients across three data frames
    common_pat = clinical$donorId[clinical$donorId %in% exp$donorId]
    common_pat = common_pat[common_pat %in% rownames(mut)]

    # Outcome vector for causal forest
    Y = clinical[common_pat, "outcome"]
    X = exp[common_pat,]

    # Count the number of patients who have a mutation in these genes
    freq_sum = sapply(mut, function(x) sum(x != 0))
    mut_selection = freq_sum > dim(mut)[1]*0.05 # select only genes with enough observations
    
    tau_p_df = data.frame()
    for(tx in which(mut_selection)) {
        W = as.numeric(mut[common_pat, tx] != 0) # treatment assignment as binary indicators
        cf = causal_forest(X, Y, W)
        cf_pred = predict(cf, estimate.variance = TRUE)
        cf_tau = compute_stats(cf_pred)
        message(paste0("The combined p values for ", colnames(mut)[tx], " is: "))
        combined_p = allmetap(cf_tau$tau.pval, method = "all") # combine tau pvalues using the metap package
        print(combined_p)
        tau_p_df = rbind(tau_p_df, c(summary(cf_tau$tau), unlist(combined_p$p) )
    }
    rownames(tau_p_df) = colnames(mut)[mut_selection]
    colnames(tau_p_df) = c(names(summary(cf_tau$tau.pval)), names(unlist(combined_p$p)))
    
    write.csv(tau_p_df, file = paste0("./exp/2020-10-07_TCGA_mixed_mut_tx_exp_covar_HTE/res/", cancer_type, "_", endpt, "_tau_pvals.csv"), row.names = TRUE)
}
