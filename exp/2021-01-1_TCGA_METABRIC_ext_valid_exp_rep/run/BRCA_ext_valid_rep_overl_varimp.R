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

##### I. PREPARE FROM RAW DATASETS

# 1. Fetch clinical data
metab_clin <- fetch_metab_clinical() # impute with simple KM method, OS in months
metab_clin <- mk_id_rownames(metab_clin)
tcga_clin <- fetch_clinical_data("BRCA", outParam = "OS", col_vec = c("bcr_patient_barcode", "type", "age_at_initial_pathologic_diagnosis", "OS", "OS.time"), outUnitDays2Month = TRUE) # tyep is still required here to subset appropriate patients from the TCGA-CDR. Also, convert to months to match with metabric data.
tcga_clin$type <-NULL
tcga_clin <-mk_id_rownames(tcga_clin)

# 2. Fetch expresion median data
metab_z <-fetch_mrna_z_score("metab", save = TRUE)
tcga_z <-fetch_mrna_z_score("tcga", save = TRUE)
tcga_z <- mk_id_rownames(tcga_z) # first colname V1
rownames(tcga_z) <- format_tcga_patient(rownames(tcga_z)) # z score data set included the sample type (i.e. "-01") even though only one type is present in the BRCA dataset
tcga_z$donorId <- NULL # still needs to be removed
metab_z <-mk_id_rownames(metab_z)


##### II. HOMOGENEIZE THE TWO DATASETS

# Subset patients with all: clinical, mutation and expression record and perform subsetting with rownames
metab_common_pat <- intersect(rownames(metab_clin), rownames(metab_z))
metab_clin <- metab_clin[metab_common_pat,]
metab_z <- metab_z[metab_common_pat,]

tcga_common_pat <- intersect(rownames(tcga_clin), rownames(tcga_z))
tcga_clin <- tcga_clin[tcga_common_pat,]
tcga_z <- tcga_z[tcga_common_pat,]

# Save list of common expression panel genes in both datasets
common_z <-intersect(colnames(metab_z), colnames(tcga_z))


 ##### III. PREPARE TX LIST AND TX DIRECTION INDICATIONS FOR BOTH DATASETS

 # 1. TX list and direction for TCGA
tcga_dea <- read.csv(paste0("./dat/tables/BRCA_DEGtable.csv"))

# Convert to SYMBOL using canonical UniProt trnascript (see cBioPortal FAQ https://docs.cbioportal.org/1.-general/faq)

# Import cBioportal conversion list
cbio_uniprot <- read.table("./raw/isoform_overrides_uniprot", sep = "\t")

# Convert ENST to ENSG
library(EnsDb.Hsapiens.v75) # Only load when script requires
edb <- EnsDb.Hsapiens.v75
esb <- AnnotationDbi::select(edb, key = cbio_uniprot$V1, keytype = "TXID", columns = c("TXID", "GENEID", "SYMBOL"))

tcga_dea_symb <- left_join(tcga_dea, esb, by = c("X" = "GENEID"), copy  = T)

# 2. same TX list and direction for METABRIC

tx_list <- intersect(tcga_dea$X)