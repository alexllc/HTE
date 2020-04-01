library(TCGAbiolinks)

cancer_types = c(
    #'ACC','BLCA','BRCA',
    'CESC','CHOL','COAD','DLBC','ESCA','GBM','HNSC','KICH', 'KIRC','KIRP','LGG','LIHC','LUAD','LUSC','MESO','OV','PAAD','PCPG','PRAD','READ', 'SKCM','STAD', 'TGCT','THCA','THYM','UCEC','UCS','UVM')

## Fixing drug names
# download.file("https://gdisc.bme.gatech.edu/Data/DrugCorrection.csv", "DrugCorrection.csv")
# drug_tbl = read.csv("DrugCorrection.csv")

# correct_drug_names <- function(df) {
#   colnames(drug_tbl)[1] <- "drug_name"
#   df <- inner_join(df, drug_tbl, by = "drug_name")
#   df$drug_name <- NULL
#   colnames(df)[colnames(df)=="Correction"] <- "drug_names"
#   return(df)
# }

# correct_drug_cat <- function(cat) {
#   targeted_rx <- read.csv(paste0(usrwd, "/cancerl/data/targeted_rx.csv"))
#   targeted_rx[,1] <- sapply(sapply(as.character(targeted_rx[,1]), function(x) strsplit(x, ' ')), function(x) x[[1]])
#   colnames(targeted_rx)[1] <- "drug_name"

# }
usrwd = "/exeh_4/alex_lau/"
setwd(paste0(usrwd, "/HTE/wd/expression_HTE"))

for (project in cancer_types) {
    print(paste0(paste0(rep('#', 20), collapse=""), "Downloading: ", project, paste0(rep('#', 20), collapse="")))
    #############################################################
    ### 1. Prepre TCGA clinical data
    #############################################################

    # Fetch survival data from GDC
    query <- GDCquery(project = paste0("TCGA-", project),
                data.category = "Clinical",
                file.type = "xml") # Obtain list from TCGAbiolinks:::getProjectSummary(<enter_project>)

    GDCdownload(query)
    # clinical_set <- c("drug", "patient")

    # for (c in clinical_set) {
    # assign(c, GDCprepare_clinic(query, clinical.info = c))
    # }

    ### Fix drug names


    #############################################################
    ### 2. Prepare 
    #############################################################

    if (!file.exists(paste0("./HTSeqData/", project, "_exp.rda")) ) {

    # Fetch expression data from GDC
    g_query <- GDCquery(project = paste0("TCGA-", project),
                    data.category = "Transcriptome Profiling",
                    legacy = F,
                    data.type = "Gene Expression Quantification",
                    workflow.type = "HTSeq - FPKM-UQ",
                    sample.type = "Primary solid Tumor")

    GDCdownload(g_query)
    expdat <- GDCprepare(query = g_query,
                        save = TRUE,
                        save.filename = paste0("./HTSeqData/",project, "_exp.rda"))
    }

}