# Functions for retreiving TCGA data in HTE analysis
## These functions requires relative path access to donwloaded files and GDCdata, therefore this script should be run under the outermost home directory.

library(TCGAbiolinks)
library(dplyr)
library(tidyr)
library(NNMIS)

impute_with_NNMIS <- function(clin_df, type = "TCGA", only_export_obj = FALSE) {
        labels = c("[Discrepancy]","[Not Applicable]","[Not Available]","[Unknown]")
        clin_df$ajcc_pathologic_tumor_stage[which(clin_df$ajcc_pathologic_tumor_stage %in% labels)] = NA
        clin_df$type = NULL
        # Conver all variables into numeric
        for (c in colnames(clin_df)) {
            if (!is.numeric(clin_df[,c]) && c != "donorId") {
                # Convert empty strings into NAs first
                which.one <- which( levels(clin_df[,c]) == "")
                levels(clin_df[,c])[which.one] <- NA
                clin_df[,c] = sapply(sapply(clin_df[,c], as.factor), as.numeric) 
                print(paste0(c, " is converted to numeric.")) 
            }
        }
        if(all(is.na(clin_df$ajcc_pathologic_tumor_stage))) {
            message("All AJCC stage entries are NA.")
            clin_df$ajcc_pathologic_tumor_stage = NULL
            clin_df$tumor_status[floor(dim(clin_df)[1]/2)] = NA # manually removing one data point or else NNMIS will not permit using this as the auxillary variable
            tcga_imp = NNMIS(clin_df$tumor_status, 
                            xa = clin_df$age_at_initial_pathologic_diagnosis, 
                            xb = clin_df$age_at_initial_pathologic_diagnosis, 
                            time = clin_df$OS.time, 
                            event = clin_df$OS, 
                            imputeCT = T, 
                            Seed = 2020, 
                            mc.cores = 60)
            tcga_imp_surv = tcga_imp$dat.T.NNMI %>% mutate(mean = rowMeans(.))
            tcga_imp_covar = as.data.frame(sapply(tcga_imp$dat.NNMI, as.numeric)) %>% mutate(mean = rowMeans(.))
            clin_df$outcome = tcga_imp_surv$mean
            clin_df$tumor_status = tcga_imp_covar$mean
        } else {
            # Removing tumor status as an additional covariate to keep as many complete cases as possible

            # Incase there are tumor types with full records we manually remove one data point or NNMIS WILL NOT RUN
            clin_df$ajcc_pathologic_tumor_stage[floor(dim(clin_df)[1]/2)] = NA

            tcga_imp = NNMIS(clin_df$ajcc_pathologic_tumor_stage, 
                            xa = clin_df$age_at_initial_pathologic_diagnosis, 
                            xb = clin_df$age_at_initial_pathologic_diagnosis, 
                            time = clin_df$OS.time, 
                            event = clin_df$OS, 
                            imputeCT = T, 
                            Seed = 2020, 
                            mc.cores = 60)
            tcga_imp_surv = tcga_imp$dat.T.NNMI %>% mutate(mean = rowMeans(.))
            tcga_imp_covar = as.data.frame(sapply(tcga_imp$dat.NNMI, as.numeric)) %>% mutate(mean = rowMeans(.))
            clin_df$outcome = tcga_imp_surv$mean
            clin_df$ajcc_pathologic_tumor_stage = tcga_imp_covar$mean
        }
    
    if (only_export_obj) {
        return(tcga_imp)
    } else {
        return(clin_df)
    }
}

#' Fetch clinical data from TCGA-CDR for building CF 
#' 
#' Fetching from (https://gdc.cancer.gov/about-data/publications/PanCan-Clinical-2018). Use this fuction in conjunction with the NNMIS_survival_imputation.R script
#'
#' @param cancer_type [string] TCGA cancer cancer_type code
#' @param col_vec [vector] of columns that need to be retreived from the TCGA-CDR, the basic requirements for NNMIS imputation is provided as default.
#' @return [dataframe] clinical info with the "donorId" column as patient code and other selected clinical columns converted to numeric
fetch_clinical_data <- function(cancer_type, col_vec =  c("bcr_patient_barcode", "type", "age_at_initial_pathologic_diagnosis",  "gender", "ajcc_pathologic_tumor_stage", "tumor_status", "OS", "OS.time", "tumor_status")) {
    if (!file.exists("./dat/TCGA-CDR-SupplementalTableS1.xlsx")) {
        download.file("https://ars.els-cdn.com/content/image/1-s2.0-S0092867418302290-mmc1.xlsx", "./dat/TCGA-CDR-SupplementalTableS1.xlsx") }

        cdr = read_excel("./dat/TCGA-CDR-SupplementalTableS1.xlsx")
        clinical_dat = dplyr::select(cdr, col_vec)

        clinical_dat[clinical_dat == "#N/A"] <- NA # sometimes missing entries are not well formatted
        clinical_dat <- subset(clinical_dat, !is.na(OS.time) & !is.na(OS) & !is.na(age_at_initial_pathologic_diagnosis))
        colnames(clinical_dat)[colnames(clinical_dat)=="bcr_patient_barcode"] = "donorId"
        clinical_dat = as.data.frame(clinical_dat)

    # specify cancer type here
    clinical_dat <- dplyr::filter(clinical_dat, type == cancer_type)
    clinical_dat = impute_with_NNMIS(clinical_dat)
    clinical_dat = dplyr::select(clinical_dat, -c(OS, OS.time))
    clinical_dat = clinical_dat[complete.cases(clinical_dat),]
    print("Processed patient dataframe: ")

    print(head(clinical_dat))

    return(clinical_dat)
}

#' Fetch mutation datat using TCGAbiolinks
#' 
#' @param cancer_type [string] TCGA cancer cancer_type code
#' @return [dataframe] with genes as columns and patient entries as rows. The first column is named "donorId" for patient id.

fetch_mut_data <- function(cancer_type) {
    m_query <- GDCquery(project = paste0("TCGA-", cancer_type),
                data.category = "Simple Nucleotide Variation",
                legacy = FALSE,
                workflow.type = "MuSE Variant Aggregation and Masking",
                data.type = "Masked Somatic Mutation"
            )
    setwd("./raw/")
    GDCdownload(m_query)
    maf = GDCprepare(m_query)
    setwd("../")

    # Cleaning and spreading MAF donwloade from TCGA
    pmaf = dplyr::filter(maf, BIOTYPE == "protein_coding")
    pmaf$SNNS = ifelse(pmaf$One_Consequence == "synonymous_variant" | is.na(pmaf$Amino_acids),"SN", "NS")
    NSpmaf = filter(pmaf, SNNS == "NS")
    pmaf = dplyr::select(pmaf, c(Hugo_Symbol, Tumor_Sample_Barcode))
    pmaf = pmaf %>% group_by(Tumor_Sample_Barcode) %>% add_count(Hugo_Symbol)
    pmaf = pmaf[!duplicated(pmaf),]
    wtcga = pmaf %>% spread(Hugo_Symbol, n)
    wtcga[is.na(wtcga)] = 0
    wtcga$Tumor_Sample_Barcode = as.character(wtcga$Tumor_Sample_Barcode)
    tmp = strsplit(wtcga$Tumor_Sample_Barcode, "-")
    tmp = unlist(lapply(tmp, function(x) paste(x[[1]], x[[2]], x[[3]], sep = "-")))
    wtcga$Tumor_Sample_Barcode <- unlist(tmp)
    wtcga = as.data.frame(wtcga)
    zero_counts = sapply(wtcga, function(x) sum(x == 0))
    colnames(wtcga)[1] = "donorId"

    return(wtcga)
}


#' Fetch RNASeq expression data from TCGAbiolinks
#' 
#' @param cancer_type [string] TCGA cancer cancer_type code
#' @param addBatch [logical] whether to return extra batch variables (they will be converted to numeric)
#' @param scale [logical] whether to scale by file maximum or not, not recommended if you are trying to plot trees
#' @return [matrix] with a "donorId" column for patient id, optionally with batch variables and ~50k genes as columns. The entries are 
#' 

fetch_exp_data <- function(cancer_type, addBatch = TRUE, scale = FLASE) {
    if (!file.exists(paste0("./dat/HTSeqData/", cancer_type, "_exp.rda")) ) {
    
    # Fetch expression data from GDC
    g_query <- GDCquery(cancer_type = paste0("TCGA-", cancer_type),
                    data.category = "Transcriptome Profiling",
                    legacy = F,
                    data.type = "Gene Expression Quantification",
                    workflow.type = "HTSeq - FPKM-UQ",
                    sample.type = "Primary Tumor")

    GDCdownload(g_query)
    expdat <- GDCprepare(query = g_query,
                        save = TRUE,
                        save.filename = paste0("./HTSeqData/",cancer_type, "_exp.rda"))
    prep <- GDCprepare(g_query) 
    } else {
        load(paste0("./dat/HTSeqData/", cancer_type, "_exp.rda"))
        prep = data
    }

    exp_matrix <- SummarizedExperiment::assay(prep, "HTSeq - FPKM-UQ")

    # Only select primary tumor samples
    # Somehow this line removed so many patients
    exp_matrix <- exp_matrix[,grep("-01.-", colnames(exp_matrix))]

    #Scale FPKM value by patients
    if (scale) exp_matrix <- exp_matrix/apply(exp_matrix,2,max)
    exp_matrix <- t(exp_matrix)

    # Address batch effects by extending tissue source, aliquot, plate and sequencing center as additional covariates
    if (add_batch) {
    exp_matrix <- cbind(separate(as.data.frame(rownames(exp_matrix)), 
                                "rownames(exp_matrix)", 
                                c(NA, "TSS", "patient", NA, "portion", "plate", "center"),  # skip var with NAs
                                sep = "-"), 
                        exp_matrix)
    }
    exp_matrix$bcr <- rownames(exp_matrix)

    # Select only one aliquot
    exp_matrix <- as.data.frame(exp_matrix %>% group_by(patient) %>% dplyr::slice(1))
    rownames(exp_matrix) <- exp_matrix$bcr


    batches = c("TSS", "patient", "portion", "plate", "center")
    for (name in batches) {
        c = grep(name, colnames(exp_matrix))
        which.one <- which( levels(exp_matrix[,c]) == "")
        levels(exp_matrix[,c])[which.one] <- NA
        print(paste0(colnames(exp_matrix)[c], " is altered"))
        exp_matrix[,c] = sapply(sapply(exp_matrix[,c], as.factor), as.numeric) 
    }

    #Format patient
    tmp = strsplit(rownames(exp_matrix), "-")
    tmp = unlist(lapply(tmp, function(x) paste(x[[1]], x[[2]], x[[3]], sep = "-")))
    rownames(exp_matrix) <- unlist(tmp)
    exp_matrix$donorId <- rownames(exp_matrix)
    exp_matrix <- dplyr::select(exp_matrix, -c(bcr, patient))

    return(exp_matrix)
}