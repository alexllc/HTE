# Perofrm external validation of HTE with a modified SHC approach using the independent datasets TCGA-BRCA and METABRIC. Both expression 
# The following parameters are compared
## tau values
## 


setwd("~/project/HTE/")

## source bin scripts
bin_ls = list.files("./bin")
for (bin in bin_ls){
    source(paste0("./bin/", bin))
}

# I. PREPARE FROM RAW DATASETS
# 1. Fetch clinical data
metab_clin <- fetch_metab_clinical() # impute with simple KM method, OS in months
metab_clin <- mk_id_rownames(metab_clin)
tcga_clin <- fetch_clinical_data("BRCA", outParam = "OS", col_vec = c("bcr_patient_barcode", "type", "age_at_initial_pathologic_diagnosis", "OS", "OS.time"), outUnitDays2Month = TRUE) # tyep is still required here to subset appropriate patients from the TCGA-CDR. Also, convert to months to match with metabric data.
tcga_clin$type <-NULL
tcga_clin <-mk_id_rownames(tcga_clin)

# 3. Fetch expresion median data
metab_z <-fetch_mrna_z_score("metab", save = TRUE)
tcga_z <-fetch_mrna_z_score("tcga", save = TRUE)
tcga_z <- mk_id_rownames(tcga_z) # first colname V1
rownames(tcga_z) <- format_tcga_patient(rownames(tcga_z)) # z score data set included the sample type (i.e. "-01") even though only one type is present in the BRCA dataset
tcga_z$donorId <- NULL # still needs to be removed
metab_z <-mk_id_rownames(metab_z)

# TODO II. HOMOGENEIZE THE TWO DATASETS
# Subset patients with all: clinical, mutation and expression record and perform subsetting with rownames