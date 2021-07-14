# Script to re-run BRCA expression HTE using imputed PFI as outcome

# Run from base directory
setwd("~/project/HTE/")

## source bin scripts
bin_ls = list.files("./bin")
for (bin in bin_ls){
    source(paste0("./bin/", bin))
}

# Clinical data
if (!file.exists(basename("TCGA_CDR_clean.csv"))) {

    download.file("https://ars.els-cdn.com/content/image/1-s2.0-S0092867418302290-mmc1.xlsx", "TCGA-CDR-SupplementalTableS1.xlsx")

    cdr = read_excel("TCGA-CDR-SupplementalTableS1.xlsx")
    clinical_dat = dplyr::select(cdr, c(bcr_patient_barcode, type, age_at_initial_pathologic_diagnosis,  gender, ajcc_pathologic_tumor_stage, PFI, PFI.time))

    #patient <- patient[patient$OS.time!=0,]
    clinical_dat[clinical_dat == "#N/A"] <- NA
    clinical_dat <- subset(clinical_dat, !is.na(PFI.time) & !is.na(PFI) & !is.na(age_at_initial_pathologic_diagnosis))
    colnames(clinical_dat)[colnames(clinical_dat)=="bcr_patient_barcode"] = "donorId"
# This renaming step is critical as the HTE main function will rely on the column named outcome to indicate Y
    #patient$gender = as.numeric(patient$gender == 'FEMALE')

    write.csv(clinical_dat, "TCGA_CDR_clean.csv", row.names = F)
} else {
    clinical_dat = read.csv("TCGA_CDR_clean.csv")
}
