# Functions for retreiving TCGA data in HTE analysis
## These functions requires relative path access to donwloaded files and GDCdata, therefore this script ***MUST*** be run under the outermost home directory.


#' Function to use an axullary variable to impute time to failture wtih a Cox PH model, must be used in conjunction with fetch_clinical_data function or the column names will not match.
#' 
#' @param clin_df [data frame] data frame with at least ("type", "age_at_initial_pathologic_diagnosis",  "gender", "ajcc_pathologic_tumor_stage", "tumor_status", outParam, paste0(outParam, ".time"), "tumor_status") as colnames, which are the defaults in the fetch_clincal_data function.
#' @param type [string] cancer type in TCGA project code
#' @param outParam [string] outcome measure to be imputed
#' @param onlyExportObj [logical] whether to return the Surv obj or the imputed dataframe. If you want both the DF and the obj, call this function twice, keep the seed the same.

impute_with_NNMIS <- function(clin_df, type = "TCGA", outParam = "OS", onlyExportObj = FALSE) {
        labels = c("[Discrepancy]","[Not Applicable]","[Not Available]","[Unknown]")
        clin_df$ajcc_pathologic_tumor_stage[which(clin_df$ajcc_pathologic_tumor_stage %in% labels)] = NA
        clin_df$type = NULL
        # Conver all variables into numeric
        clin_df = convert_col_to_numeric(clin_df)
        if(all(is.na(clin_df$ajcc_pathologic_tumor_stage))) {
            message("All AJCC stage entries are NA.")
            clin_df$ajcc_pathologic_tumor_stage = NULL
            clin_df$tumor_status[floor(dim(clin_df)[1]/2)] = NA # manually removing one data point or else NNMIS will not permit using this as the auxillary variable
            tcga_imp = NNMIS(clin_df$tumor_status, 
                            xa = clin_df$age_at_initial_pathologic_diagnosis, 
                            xb = clin_df$age_at_initial_pathologic_diagnosis, 
                            time = clin_df$out_param_time, 
                            event = clin_df$out_param, 
                            imputeCT = T, 
                            Seed = 2020, 
                            mc.cores = 60)
            tcga_imp_surv = tcga_imp$dat.T.NNMI %>% mutate(mean = rowMeans(.))
            tcga_imp_covar = as.data.frame(sapply(tcga_imp$dat.NNMI, as.numeric)) %>% mutate(mean = rowMeans(.))
            clin_df$outcome = tcga_imp_surv$mean
            clin_df$tumor_status = tcga_imp_covar$mean
        } else {

            clin_df$ajcc_pathologic_tumor_stage[floor(dim(clin_df)[1]/2)] = NA # Incase there are tumor types with full records we manually remove one data point or NNMIS WILL NOT RUN
            clin_df$tumor_status = NULL # Removing tumor status as an additional covariate to keep as many complete cases as possible
            tcga_imp = NNMIS(clin_df$ajcc_pathologic_tumor_stage, 
                            xa = clin_df$age_at_initial_pathologic_diagnosis, 
                            xb = clin_df$age_at_initial_pathologic_diagnosis, 
                            time = clin_df$out_param_time, 
                            event = clin_df$out_param, 
                            imputeCT = T, 
                            Seed = 2020, 
                            mc.cores = 60)
            tcga_imp_surv = tcga_imp$dat.T.NNMI %>% mutate(mean = rowMeans(.))
            tcga_imp_covar = as.data.frame(sapply(tcga_imp$dat.NNMI, as.numeric)) %>% mutate(mean = rowMeans(.))
            clin_df$outcome = tcga_imp_surv$mean
            clin_df$ajcc_pathologic_tumor_stage = tcga_imp_covar$mean
        }
    
    if (onlyExportObj) {
        return(tcga_imp)
    } else {
        return(clin_df)
    }
}

#' Fetch clinical data from TCGA-CDR for building CF 
#' 
#' Fetching from (https://gdc.cancer.gov/about-data/publications/PanCan-Clinical-2018). Use this fuction in conjunction with the NNMIS_survival_imputation.R script
#'
#' @param cancer_type 
#' 
#' [string] TCGA cancer cancer_type code
#' @param outParam [string] outcome parameter to be examined, default is OS, but you can select from OS, PFI, DFI or DSS.
#' @param outUnitDays2Month whether to convert outcome unit from days to months
#' @param imputeMethod [string] ways to impute outcome length missing due to censoring issues, default is a simple KM plot.
#' @param onlyCompleteCases [bool] whether to return complete cases only, not recommended for prognostic models in general
#' @param col_vec [vector] of columns that need to be retreived from the TCGA-CDR, the basic requirements for NNMIS imputation is provided as default.
#' @return [dataframe] clinical info with the "donorId" column as patient code and other selected clinical columns converted to numeric
fetch_clinical_data <- function(cancer_type, outParam = "OS", outUnitDays2Month = FALSE, imputeMethod = "simple", onlyCompleteCases = FALSE, col_vec =  c("bcr_patient_barcode", "type", "age_at_initial_pathologic_diagnosis",  "gender", "ajcc_pathologic_tumor_stage", "tumor_status", outParam, paste0(outParam, ".time"))) {
    if (!file.exists("./dat/TCGA-CDR-SupplementalTableS1.xlsx")) {
        download.file("https://ars.els-cdn.com/content/image/1-s2.0-S0092867418302290-mmc1.xlsx", "./dat/TCGA-CDR-SupplementalTableS1.xlsx") }

        cdr = read_excel("./dat/TCGA-CDR-SupplementalTableS1.xlsx")
        clinical_dat = dplyr::select(cdr, col_vec)

        clinical_dat[clinical_dat == "#N/A"] <- NA # sometimes missing entries are not well formatted
        clinical_dat <- subset(clinical_dat, !is.na(get(paste0(outParam, ".time"))) & !is.na(get(outParam)) & !is.na(age_at_initial_pathologic_diagnosis))
        colnames(clinical_dat)[colnames(clinical_dat)=="bcr_patient_barcode"] = "donorId"
        colnames(clinical_dat)[colnames(clinical_dat) == outParam | colnames(clinical_dat) == paste0(outParam, ".time")] = c("out_param", "out_param_time")
        clinical_dat = as.data.frame(clinical_dat)

    # specify cancer type here
    clinical_dat <- dplyr::filter(clinical_dat, type == cancer_type)
    
    if (imputeMethod == "simple") {
        imput_surv = exp(impute.survival(clinical_dat$out_param_time, clinical_dat$out_param))
        if (outUnitDays2Month) imput_surv = imput_surv / 30.417 # reverse log and convert days to months for better interpretation. Converting ratio 30.417 is an approximate value used by Google.
        clinical_dat$outcome = imput_surv
        # Variable type conversion is not intrinsically implemented in the impute.survival function, so it must implemented here.
        clinical_dat = convert_col_to_numeric(clinical_dat)
    } else {
        clinical_dat = impute_with_NNMIS(clinical_dat)
    }
    clinical_dat = dplyr::select(clinical_dat, -c(out_param, out_param_time))
    if (onlyCompleteCases) clinical_dat = clinical_dat[complete.cases(clinical_dat),]
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
    GDCdownload(m_query)  # can't explicitly set directory, requires setwd
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
#' @param addBatch [bool] whether to return extra batch variables (they will be converted to numeric)
#' @param numericBatch [bool] convert batch factor levels into numeric. Useful for using these variables as 
#' @param scale [bool] whether to scale by file maximum or not, not recommended if you are trying to plot trees
#' @param primaryTumorOnly [bool] option to select whether to output a dataframe with other non-primary tissue expression entries. By default you should not format the tissue ID to avoid duplicates.
#' @return [matrix] with a "donorId" column for patient id, optionally with batch variables and ~50k genes as columns. The entries are 
#' 

fetch_exp_data <- function(cancer_type, addBatch = TRUE, numericBatch = TRUE, scale = FALSE, primaryTumorOnly = FALSE, formatPatient = FALSE) {
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
    if (primaryTumorOnly) exp_matrix <- exp_matrix[,grep("-01.-", colnames(exp_matrix))]

    #Scale FPKM value by patients
    if (scale) exp_matrix <- exp_matrix/apply(exp_matrix,2,max)
    exp_matrix <- t(exp_matrix)

    # Address batch effects by extending tissue source, aliquot, plate and sequencing center as additional covariates
    if (addBatch) {
        exp_matrix = cbind(separate(as.data.frame(rownames(exp_matrix)), 
                                    "rownames(exp_matrix)", 
                                    c(NA, "TSS", "patient", NA, "portion", "plate", "center"),  # skip var with NAs
                                    sep = "-"), 
                            exp_matrix)
    }

    # Select only one aliquot
    message(paste0("Expression matrix size went from :", dim(exp_matrix)[1], " * ", dim(exp_matrix)[2]))
    prev_sample = rownames(exp_matrix)
    if (primaryTumorOnly) exp_matrix <- as.data.frame(exp_matrix %>% rownames_to_column("donorId") %>%group_by(patient) %>% dplyr::slice(1) %>% column_to_rownames("donorId"))
    message(paste0("to :", dim(exp_matrix)[1], " * ", dim(exp_matrix)[2]))
    # Report duplicate
    miss_sample = prev_sample[!(prev_sample %in% rownames(exp_matrix))]
    for (sample in miss_sample){
        matches = strsplit(sample, "-")[[1]][3]
        prev_match = prev_sample[grep(matches, prev_sample)]
        selected_match = rownames(exp_matrix)[grep(matches, rownames(exp_matrix))]
        prev_match = prev_match 
    }

    batches = c("TSS", "patient", "portion", "plate", "center")

    # Convert all batch covaraites into numeric
    if (addBatch & numericBatch) {
        for (name in batches) {
            c = grep(name, colnames(exp_matrix))
            which.one <- which( levels(exp_matrix[,c]) == "")
            levels(exp_matrix[,c])[which.one] <- NA
            print(paste0(colnames(exp_matrix)[c], " is converted to numeric"))
            exp_matrix[,c] = sapply(sapply(exp_matrix[,c], as.factor), as.numeric) 
        } 
    } else if(!addBatch & numericBatch) {
            warning("No batch varaible added, not converted to numeric.")
        }

    if(primaryTumorOnly & formatPatient) {
        rownames(exp_matrix) = format_tcga_patient(rownames(exp_matrix))
    } else if (!primaryTumorOnly & formatPatient) {
       warning("Formatting tumor barcode without first selecting for primary tumor will result in duplicated entries, not formatting patient barcode.")
    }

    return(exp_matrix)
}

mk_id_rownames <- function(df) {
    rownames(df) = df$donorId
    df$donorId = NULL
    return(df)
}

convert_col_to_numeric <- function(clin_df) {
    for (c in colnames(clin_df)) {
        if (!is.numeric(clin_df[,c]) && c != "donorId") {
            # Convert empty strings into NAs first
            which.one <- which( levels(clin_df[,c]) == "")
            levels(clin_df[,c])[which.one] <- NA
            clin_df[,c] = sapply(sapply(clin_df[,c], as.factor), as.numeric) 
            print(paste0(c, " is converted to numeric.")) 
        }
    }
    return(clin_df)
}


#' Only keep the first three TCGA barcode item: project, TSS and participants. This is useful for merging with clinical data frames, but not for tissue-based data.
#' 
#' @param pat_ls [string vector] list with strings of TCGA barcode

format_tcga_patient <- function(pat_ls) {
    tmp = strsplit(pat_ls, "-")
    tmp = unlist(lapply(tmp, function(x) paste(x[[1]], x[[2]], x[[3]], sep = "-")))
    return(unlist(tmp))
}


#' Function to filter technical replicates in TCGA samples
#' Doc: exp/2020-10-07_TCGA_mixed_mut_tx_exp_covar_HTE/doc/2020-11-12_technical_replicates_in_TCGA.md
#' 
#' @source \url{http://gdac.broadinstitute.org/runs/stddata__2014_01_15/samples_report/READ_Replicate_Samples.html}
#' @param bcr list of barcodes
#' 
#' @return subset of the list of barcodes
#' 
filter_replicate_samples <- function(bcr) {
    if ( all(grepl(".{19}[RHT]", bcr)) ) { 
        type <- "RNA"
    } else if (all (grepl(".{19}[DGWX]"))) {
       type <- "DNA"
    } else {
        stop("Mixing RNA and DNA samples not allowed.")
    }

    old_len <- length(bcr)

    bcr_df <- cbind(separate(as.data.frame(bcr),
                                    "bcr",
                                    c("project", "TSS", "patient", "sample,vial", "portion,analyte", "plate", "center"),  # skip var with NAs
                                    sep = "-"),
                                    bcr)
    bcr_df <- separate(bcr_df, col = "sample,vial", into = c("sample", "vial"), sep = 2)
    bcr_df <- separate(bcr_df, col = "portion,analyte", into = c("portion", "analyte"), sep = 2)
    bcr_df <- bcr_df %>% arrange(analyte, desc(plate), desc(portion)) %>% group_by(TSS, patient, sample) %>% slice(1)

    new_len <- length(bcr_df$bcr)

    if (new_len == old_len) {
        print("No sample needs to be removed.")
        return(bcr_df$bcr)
    } else {
        out_tbl <- data.frame()
        dup_samples <- unique(substr(bcr[!(bcr %in% bcr_df$bcr)], 1, 15))
        for (dupbcr in dup_samples) {
            kept <- as.character(bcr_df[grep(dupbcr, bcr_df$bcr), 10])
            removed <- bcr[grep(dupbcr, bcr)]
            removed <- paste(removed[!(removed %in% kept)], collapse = ", ")
            out_tbl <- rbind(out_tbl, c(kept, removed))
        }
        colnames(out_tbl) <- c("chosen", "removed")
        print(paste0("removed the following samples: ", out_tbl$removed))
        # out_tbl
        return(bcr_df$bcr)
    }
}