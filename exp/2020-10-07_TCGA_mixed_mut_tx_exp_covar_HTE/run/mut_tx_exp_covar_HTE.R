# Script for running HTE analysis using mutation as treatment variables and expression as covaraites 

## Run from base directory
setwd("~/project/HTE/")

## source bin scripts
bin_ls = list.files("./bin")
for (bin in bin_ls){
    source(paste0("./bin/", bin))
}

## Set cacner type for this run
cancer_type = "BRCA"
endpt = "OS"

## Prepare mutation dataframe
# Loading the previously wide-formatted MAF is preferred as re-formatting takes too long
if(!file.exists(paste0("./dat/mutation_frequencies/MAF_wide_", cancer_type, ".csv.gz"))) {
    mut = as.data.frame(fetch_mut_data(cancer_type)) # will take a long time, you'd want to save it
    write.csv(mut, file = gzfile(paste0("./dat/mutation_frequencies/MAF_wide_", cancer_type, ".csv.gz")), row.names = FALSE)
} else {
    mut = as.data.frame(fread(paste0("./dat/mutation_frequencies/MAF_wide_", cancer_type, ".csv.gz")))
}
mut = mk_id_rownames(mut)

## Prepare clinical dataframe
clinical = fetch_clinical_data(cancer_type, outParam = endpt, imputeMethod = "simple", outUnitDays2Month = TRUE)
clinical = mk_id_rownames(clinical)

## Prepare expression dataframe
exp = fetch_exp_data(cancer_type, scale = FALSE, primaryTumorOnly = TRUE, formatPatient = TRUE)
exp = mk_id_rownames(exp)

## Subset patients and settings for HTE
# Check the list of common patients across three data frames
common_pat = rownames(clinical)[rownames(clinical) %in% rownames(exp)]
common_pat = common_pat[common_pat %in% rownames(mut)]

# Outcome vector for causal forest
Y = clinical[common_pat, "outcome"]
X = exp[common_pat,]
# Count the number of patients who have a mutation in these genes
freq_sum = sapply(mut, function(x) sum(x != 0))
mut_selection = freq_sum > dim(mut)[1]*0.05 # select only genes with enough observations


## Run and save HTE