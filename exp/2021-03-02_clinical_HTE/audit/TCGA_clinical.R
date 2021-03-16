library(TCGAbiolinks)
query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab")
GDCdownload(query)
clinical.BCRtab.all <- GDCprepare(query)


names(clinical.BCRtab.all)

# [1] "clinical_follow_up_v2.1_brca"     "clinical_follow_up_v4.0_nte_brca"
# [3] "clinical_radiation_brca"          "clinical_omf_v4.0_brca"          
# [5] "clinical_patient_brca"            "clinical_follow_up_v4.0_brca"    
# [7] "clinical_drug_brca"               "clinical_nte_brca"               
# [9] "clinical_follow_up_v1.5_brca"  

# nte = new tumor event
# omf = other malignancy form

for (tbl in names(clinical.BCRtab.all)) {
    print(paste(rep("=", 50), collapse = ""))
    print(tbl)
    print(paste(rep("=", 50), collapse = ""))
    # print(sort(colnames(clinical.BCRtab.all[[tbl]])))
    
    # assign(tbl, clinical.BCRtab.all[[tbl]])

    print(as.data.frame(filter(clinical.BCRtab.all[[tbl]], bcr_patient_barcode == missing_OS)))
}

nte <- clinical.BCRtab.all$clinical_nte_brca
fu1.5 <- clinical.BCRtab.all$clinical_follow_up_v1.5_brca
fu2.1 <- clinical.BCRtab.all$clinical_follow_up_v2.1_brca


missing_OS <- treatment_W[which(!(names(treatment_W) %in% clin_covar$donorId))]

TCGA-OL-A66H

# Check blatant errors in survival days of clinical indexed data
days_col <- colnames(clinical)[grep("days_to", colnames(clinical))]
days_col <- colnames(clin_BCR)[grep("days_to", colnames(clin_BCR))]

for (days in days_col){
    print(paste(rep("=", 50), collapse = ""))
    print(days)
    print(paste(rep("=", 50), collapse = ""))

    print(summary(as.numeric(clin_BCR[[days]])))
}


## Missing stages
no_stage <- clin_covar$bcr_patient_barcode[which(is.na(clin_covar$ajcc_pathologic_tumor_stage))]

find_clincal <- function(bcr = NULL) {
    for (tbl in names(clinical.BCRtab.all)) {
        if (bcr %in% clinical.BCRtab.all[tbl][["bcr_patient_barcode"]]) {
            print(paste(rep("=", 50), collapse = ""))
            print(filter(tbl, bcr_patient_barcode == bcr))
            print(paste(rep("=", 50), collapse = ""))
        } else {
            message("not here.")
        }
    }
}