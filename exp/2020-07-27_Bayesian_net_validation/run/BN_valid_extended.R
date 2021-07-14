# revised unctional validation of genes identified by Bayesian network

file_path <- "/home/alex/project/HTE/exp/2020-07-27_Bayesian_net_validation/raw/Gene_Lists"
gene_files <- list.files(file_path)

conditions <- c("Breast", "Prostate")

for (dx in conditions) {
    dx_files <- gene_files[grep(dx, gene_files)]
    bg_files <- gene_files[grep("Background", dx_files)]

# enrichement analysis

# survival relevance

# Significance in GRF