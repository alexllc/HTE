#' Script to perform external validation with METABRIC 
# Run from base directory
setwd("~/project/HTE/")

# source bin scripts
bin_ls = list.files("./bin")
for (bin in bin_ls){
    source(paste0("./bin/", bin))
}

# Parameters for this run
is_tuned = 
is_binary = 



# Download METABRIC and TCGA files from cBioportal if not already done so, these two conditions prevents script progression if raw files not available, you will need to extract 
if( !file.exists("./raw/cbioportal_METABRIC/")) {
    download.file("http://download.cbioportal.org/brca_metabric.tar.gz", "./raw/")
    untar("./raw/brca_metabric.tar.gz", exdir = "./dat/METABRIC-TCGA_external_validation/")
}

if (!file.exists("./dat/METABRIC-TCGA_external_validation/")){
    download.file("http://download.cbioportal.org/brca_tcga_pan_can_atlas_2018.tar.gz", "./raw/")
    untar("./raw/‘brca_tcga_pan_can_atlas_2018.tar.gz.1’", exdir = "./dat/METABRIC-TCGA_external_validation/")
}

# Fetch clinical data
metab_clin = fetch_metab_clinical()
metab_clin = mk_id_rownames(metab_clin)
tcga_clin = fetch_clinical_data("BRCA", outParam = "OS", col_vec = c("bcr_patient_barcode", "type", "age_at_initial_pathologic_diagnosis", "OS", "OS.time")) # tyep is still required here to subset appropriate patients from the TCGA-CDR
tcga_clin$type = NULL
tcga_clin = mk_id_rownames(tcga_clin)

# Fetch mutation data
metab_mut = fetch_metab_mut()
metab_mut = mk_id_rownames(metab_mut)

tcga_mut = as.data.frame(fread("./dat/mutation_frequencies/MAF_wide_BRCA.csv.gz"))
tcga_mut = mk_id_rownames(tcga_mut)

common_mut = intersect(colnames(metab_mut), colnames(tcga_mut))


# Fetch expresion median data
metab_z = fetch_mrna_z_score("metab", save = TRUE)
metab_z = mk_id_rownames(metab_z)

tcga_z = fetch_mrna_z_score("tcga", save = TRUE)
metab_z = mk_id_rownames(metab_z)

tcga_z$donorId = format_tcga_patient(tcga_z$donorId)
common_z = intersect(colnames(metab_z), colnames(tcga_z))
common_z = common_z[-which(common_z == "donorId")]

tcga_common_pat <- 

X <- tcga_z[]

for(tx in common_mut) {
    

    # Build individual forest
    tcga_cf = causal_forest(
                            X <- tcga_z,
                            Y <- ,
                            W,
                            Y.hat = NULL,
                            W.hat = NULL,
                            num.trees = 2000,
                            sample.weights = NULL,
                            clusters = NULL,
                            equalize.cluster.weights = FALSE,
                            sample.fraction = 0.5,
                            mtry = min(ceiling(sqrt(ncol(X)) + 20), ncol(X)),
                            min.node.size = 5,
                            honesty = TRUE,
                            honesty.fraction = 0.5,
                            honesty.prune.leaves = TRUE,
                            alpha = 0.05,
                            imbalance.penalty = 0,
                            stabilize.splits = TRUE,
                            ci.group.size = 2,
                            tune.parameters = "none",
                            tune.num.trees = 200,
                            tune.num.reps = 50,
                            tune.num.draws = 1000,
                            compute.oob.predictions = TRUE,
                            orthog.boosting = FALSE,
                            num.threads = NULL,
                            seed = runif(1, 0, .Machine$integer.max)
                        )
    # 
}