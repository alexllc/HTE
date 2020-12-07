# Script for running HTE analysis using mutation as treatment variables and expression as covaraites 

## Run from base directory
setwd("~/project/HTE/")

## source bin scripts
bin_ls = list.files("./bin")
for (bin in bin_ls){
    source(paste0("./bin/", bin))
}

## Set parameters for this run
cancer_type = "BRCA"
endpt = "OS"
output_file = "./exp/2020-10-07_TCGA_mixed_mut_tx_exp_covar_HTE/res/2020-11-30_perm_all/"
mut_thres = 0.05

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
clinical = fetch_clinical_data(cancer_type, outParam = endpt, imputeMethod = "simple", outUnitDays2Month = TRUE, discard = c("type", "tumor_status"))
clinical
clinical = mk_id_rownames(clinical)

## Prepare expression dataframe
exp = fetch_exp_data(cancer_type, scale = FALSE, primaryTumorOnly = TRUE, formatPatient = TRUE)

## Subset patients and settings for HTE
# Check the list of common patients across three data frames
common_pat = rownames(clinical)[rownames(clinical) %in% rownames(exp)]
common_pat = common_pat[common_pat %in% rownames(mut)]
message(paste0("Number of patients with expression, mutation and clinical entries is: ", length(common_pat)))

# Outcome vector for causal forest
Y = clinical[common_pat, "outcome"]
X = exp[common_pat,]
# Count the number of patients who have a mutation in these genes
freq_sum = sapply(mut, function(x) sum(x != 0))
mut_selection = freq_sum > dim(mut)[1]*mut_thres # select only genes with enough observations
tx_list = colnames(mut)[mut_selection]
X = inner_join(rownames_to_column(mut[common_pat, tx_list]), rownames_to_column(X), by = "rowname")

whole_dat = inner_join(rownames_to_column(clinical), X, by = "rowname", suffix = c("", ""))
X$rowname = NULL
colnames(whole_dat)[colnames(whole_dat) == "rowname"] = "donorId"

# Assign treatment group
W = create_tx_matrix(txVector = tx_list, 
                    binaryVector = rep(TRUE, length(tx_list)), 
                    cutoffThreshDf = data.frame(
                                        dirct = rep(">", length(tx_list)), 
                                        thresh = rep(0, length(tx_list))
                                        ), 
                    covarMat = X)

## Run and save HTE
obsNumber <- dim(X)[1]
trainId <- sample(1: obsNumber, floor(obsNumber/2), replace = FALSE)
registerDoParallel(10)

result <- run.hte(covar_mat = X, 
                tx_vector = tx_list, 
                whole_dataset = whole_dat, 
                project = cancer_type, 
                W_matrix = W, 
                trainId = trainId, 
                seed = 111, 
                is_binary = T, 
                is_save = T, 
                save_split = T, 
                is_tuned = F, 
                thres = 0.75, 
                n_core = 6, 
                output_directory = output_file, 
                perm_all = TRUE) # pre-filtered thres, should not have an effect here

write.csv(result[[1]], paste0(output_file, cancer_type, '_mixed_HTE_correlation_test_result.csv'), quote = F, row.names = F)
write.csv(result[[2]], paste0(output_file, cancer_type, '_mixed_HTE_calibration_result.csv'), quote = F, row.names = F)
write.csv(result[[3]], paste0(output_file, cancer_type, '_mixed_HTE_median_t_test_result.csv'), quote = F, row.names = F)
write.csv(result[[4]], paste0(output_file, cancer_type, '_mixed_HTE_permutate_testing_result.csv'), quote = F, row.names = F)
