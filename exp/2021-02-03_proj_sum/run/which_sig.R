library(dplyr)

cancerList <- c("BLCA", "BRCA", "COAD", "GBM", "HNSC", "LUAD", "LUSC", "KIRC", "PRAD", "STAD", "THCA", "UCEC")

sum_extract <- NULL

for (cancer in cancerList) {
    # sum_file <- read.csv(paste0("2020-06-20_TCGA_pancancer_DEA_indicated_tx_dirct_.25_gp_as_tx_exp_HTE_", cancer, "_summary.csv"))
    sum_file <- read.csv(paste0("2020-08-31_TCGA_pancancer_pathway_PROPS_HTE_", cancer, "_summary.csv"))
    # sum_file <- filter(sum_file, diff.pred.pval < 0.05 & mean.pred.pval < 0.05)
    sum_file <- filter(sum_file, diff.pred.pval < 0.05)
    tissue <- rep(cancer, dim(sum_file)[1])
    sum_file <- cbind(tissue, sum_file)
    sum_extract <- rbind(sum_extract, sum_file)
    colnames(sum_extract) <- colnames(sum_file)
}

write.csv(sum_extract, "pathway_unfiltered_summary.csv", row.names = FALSE)

