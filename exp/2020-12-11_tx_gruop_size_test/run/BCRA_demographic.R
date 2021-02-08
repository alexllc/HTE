setwd("~/project/HTE/")

## source bin scripts
bin_ls = list.files("./bin")
for (bin in bin_ls){
    source(paste0("./bin/", bin))
}

# Import TCGA clinical data
tcga_clin <- fread("./exp/2020-12-11_tx_gruop_size_test/dat/735bc5ff-86d1-421a-8693-6e6f92055563/nationwidechildrens.org_clinical_patient_brca.txt")
legacy_col <- c("bcr_patient_barcode", "birth_days_to", "gender", "menopause_status", "er_status_by_ihc", "her2_status_by_ihc")
tcga_clin <- dplyr::select(tcga_clin, all_of(legacy_col))
tcga_clin <- tcga_clin[-c(1:2),]
tcga_clin$age <- as.numeric(tcga_clin$birth_days_to) * -0.00273973
table(cut(tcga_clin$age, breaks=c(20, 30, 40, 50, 60, 70, 80, 90, 100), right = FALSE))

# Import METABRIC clinical data
metab_clin <- fread("./raw/cbioportal_METABRIC/data_clinical_patient.txt", skip = 4)
metab_col <- c("PATIENT_ID", "AGE_AT_DIAGNOSIS", "ER_IHC", "INFERRED_MENOPAUSAL_STATE")
metab_clin <- dplyr::select(metab_clin, all_of(metab_col))