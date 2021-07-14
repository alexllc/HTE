# working script to corrrect drugs used in all 33 types of cancer in TCGA

library(TCGAbiolinks)

# import all drug names registered on TCGAbiolinks
projects <- getGDCprojects()
projects <- projects$id[grep("TCGA", projects$id)]

query <- GDCquery(project = projects, 
                data.category = "Clinical",
                data.type = "Clinical Supplement", 
                data.format = "BCR Biotab")
GDCdownload(query)
clinical.BCRtab.all <- GDCprepare(query)
names(clinical.BCRtab.all)
drug <- clinical.BCRtab.all[[paste0("clinical_drug_", tolower(cancerType))]]