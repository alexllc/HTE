#' Script to perform external validation with METABRIC 
# Run from base directory
setwd("~/project/HTE/")

# source bin scripts
bin_ls = list.files("./bin")
for (bin in bin_ls){
    source(paste0("./bin/", bin))
}

# Download METABRIC and TCGA files from cBioportal if not already done so, these two conditions prevents script progression if raw files not available, you will need to extract 
if( !file.exists("./raw/cbioportal_METABRIC/")) {
    download.file("http://download.cbioportal.org/brca_metabric.tar.gz", "~/raw/")
    untar("./raw/brca_metabric.tar.gz", exdir = )
}

if (!file.exists("./dat/METABRIC-TCGA_external_validation/")){
    download.file()
}

# Fetch clinical data
metab_clin = fetch_metab_clinical()
tcga_clin = fetch_clinical_data("BRCA", outParam = "OS", col_vec = c("bcr_patient_barcode", "type", "age_at_initial_pathologic_diagnosis", "OS", "OS.time")) # tyep is still required here to subset appropriate patients from the TCGA-CDR
tcga_clin$type = NULL

# Fetch mutation data
metab_mut = fetch_metab_mut()
tcga_mut = as.data.frame(fread("./dat/mutation_frequencies/MAF_wide_BRCA.csv.gz"))
common_mut = intersect(colnames(metab_mut), colnames(tcga_mut))


# Fetch expresion median data
metab_z = fetch_mrna_z_score("metab", save = TRUE)
tcga_z = fetch_mrna_z_score("tcga", save = TRUE)

# Common entries in all three dataframes


# Perform 