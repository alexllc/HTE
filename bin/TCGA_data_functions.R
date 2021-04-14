# Functions for retreiving TCGA data in HTE analysis
## These functions requires relative path access to donwloaded files and GDCdata, therefore this script ***MUST*** be run under the outermost home directory.

# TODO update survival imputation methods in accordance to NICE TA TSD14 guidelines (http://nicedsu.org.uk/technical-support-documents/survival-analysis-tsd/)

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
#' @param discard clinical columns to get rid of
#' @return [dataframe] clinical info with the "donorId" column as patient code and other selected clinical columns converted to numeric
fetch_clinical_data <- function(cancer_type, 
                                outParam = "OS", 
                                outUnitDays2Month = FALSE, 
                                imputeMethod = "simple", 
                                onlyCompleteCases = FALSE, 
                                col_vec =  c("bcr_patient_barcode", "type", "age_at_initial_pathologic_diagnosis",  "gender", "ajcc_pathologic_tumor_stage", "tumor_status", outParam, paste0(outParam, ".time")), 
                                discard = c("type")) {
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
    clinical_dat = dplyr::select(clinical_dat, -c(out_param, out_param_time, discard))
    if (onlyCompleteCases) clinical_dat = clinical_dat[complete.cases(clinical_dat),]
    print("Processed patient dataframe: ")

    print(head(clinical_dat))

    return(clinical_dat)
}

# TODO scavenge clinical information from up-to-date GDC data
#' Fetch up to date clinical data
#' 
#' Although clinical outcome in the TCGA-CDR was more standardized calculation, it's more up-to-date and we'd sometimes want to have more cases rather than being accurate alone.
#' This function is written for OS as endpoint only, other endpoins have not been implemented yet. 
#' "OS is the period from the date of diagnosis until the date of death from any cause. "

find_clinical_data <- function(cancerType, col_vec =  c("bcr_patient_barcode", "type", "age_at_initial_pathologic_diagnosis",  "gender", "ajcc_pathologic_tumor_stage")) {
    
    # Check if TCGA-CDR is missing some patients
    cdr <- read_excel("./dat/TCGA-CDR-SupplementalTableS1.xlsx")
    cdr <- filter(cdr, type == cancerType)
    
    # Check if newer patients have been added to the clinical indexed files
    clinical <- GDCquery_clinic(project = paste0("TCGA-", cancerType), type = "clinical")
    if (!all(clinical$bcr_patient_barcode %in% cdr$bcr_patient_barcode)) {
        # CDR is not updated, some patients are missing from CDR
        new_patients <- clinical[!(clinical$bcr_patient_barcode %in% cdr$bcr_patient_barcode),]
    }
}

fetch_BCR_clinical <- function(cancerType) {

    # clin data via BCR biotab files
    query <- GDCquery(project = paste0("TCGA-", cancerType), 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab")
    GDCdownload(query)
    clinical.BCRtab.all <- GDCprepare(query)
    clin_BCR <- clinical.BCRtab.all$clinical_patient_brca
    clin_covar <- dplyr::select(clin_BCR, all_of(c("bcr_patient_barcode", "age_at_diagnosis", "ajcc_pathologic_tumor_stage", "gender", "vital_status", "last_contact_days_to", "death_days_to", "days_to_initial_pathologic_diagnosis")))
    clin_covar <- clin_covar[-c(1:2),] # remove col IDs and descriptions
    clin_covar[clin_covar == "[Not Available]"] <- NA # remove string NA identifiers
    numeric_cols <- grep("age_|days_", colnames(clin_covar))
    clin_covar[numeric_cols] <- lapply(clin_covar[numeric_cols], as.numeric)
    # convert to numeric

    # Patients are either: ALIVE, DEAD or NA in this data. OS lengths of alive patients were indicated by last_contact_days_to, while dead patients were indicated by death_days_to. Once you have one entry, the other one should be NA, so to get the "OS length" you will need to coalesce these two cols
    clin_covar <- clin_covar %>% mutate(OS_time = coalesce(last_contact_days_to, death_days_to))

    # set negative OS times to 0s, these patients were contacted before diagnosis was made, and perhaps had been lost to follow-up
    clin_covar$OS_time[clin_covar$OS_time < 0] <- 0
    clin_covar <- dplyr::select(clin_covar, -c("last_contact_days_to", "death_days_to", "days_to_initial_pathologic_diagnosis"))

    # Convert categorical variables into numeric for imputing missing stages
    clin_covar <- convert_col_to_numeric(as.data.frame(clin_covar), id = "bcr_patient_barcode")

    # impute missing stages with kNN
    tmp = kNN(as.matrix(dplyr::select(clin_covar, -c("bcr_patient_barcode", "OS_time"))))

    # patient missing survival data
    missing_surv <- clin_covar[is.na(clin_covar$days_to_last_follow_up),]

    
    # Convert "days_to_last_follow_up" negative values to 0s, these patients were recorded before a proper diagnosis was made


    # impute missing stages with 
}

#' Fetch clinical data from the indexed clinial files from TCGAbiolinks. Refer to the `fetch_BCR_clinical` function if you wish to source your clinical data from the BCRbio tabs instead. The distinction between these files are docuemnted in the "Useful Information" box in the TCGAbiolinks clinical documentation. By default, this function will try to compute missing tumor stages based on T/N/M indicators using the `find_ajcc_stage` function. However, since GBM patients do not use such staging system, the indexed clinical dataframe output will omitt the `ajcc_pathologic_stage` column.
#' 
#' @param cancerType [string] TCGA project code
#' 
#' @source \url{https://www.bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/clinical.html}
#' @return [dataframe] with the following columns: "bcr_patient_barcode", "age_at_index", "ajcc_pathologic_stage", "OS_time"
#' 
fetch_indexed_clinical <- function(cancerType) {
    clinical <- GDCquery_clinic(project = paste0("TCGA-", cancerType), type = "clinical")
    # days_to_last_known_disease_status is all NAs

    clinical_sub <- try(dplyr::select(clinical, all_of(c("bcr_patient_barcode", "age_at_index", "ajcc_pathologic_stage", "ajcc_pathologic_t", "ajcc_pathologic_n", "ajcc_pathologic_m", "days_to_diagnosis", "days_to_last_follow_up", "days_to_death", "vital_status"))))
    # PRAD is known to have a missing ajcc_pathological_m column, in this case we will use the ajcc_clinical_m column as a replacement
    if (class(clinical_sub) == "try-error") {
        message("Some AJCC pathologic M indicator missing, using AJCC clinical M instead.")
        clinical_sub <- dplyr::select(clinical, all_of(c("bcr_patient_barcode", "age_at_index", "ajcc_pathologic_t", "ajcc_pathologic_n", "ajcc_clinical_m", "days_to_diagnosis", "days_to_last_follow_up", "days_to_death", "vital_status")))
        colnames(clinical_sub)[colnames(clinical_sub) == "ajcc_clinical_m"]  <- "ajcc_pathologic_m"
        clinical_sub$ajcc_pathologic_stage <- NA
    }
    clinical_sub <- clinical_sub %>% mutate(OS_time = coalesce(days_to_last_follow_up, days_to_death))
    if (cancerType != "GBM") {
        missing_stages <- which(is.na(clinical_sub$ajcc_pathologic_stage))

        # Omit samples with more than 1 indicators
        tnm_sub <- dplyr::select(clinical_sub, all_of(c("ajcc_pathologic_t", "ajcc_pathologic_n", "ajcc_pathologic_m")))[missing_stages,]
        clinical_sub <- clinical_sub[-as.numeric(names(which(apply(tnm_sub, 1, function(x) all(is.na(x)))))),]

        # See if you can re-assign patients with missing stages using the TNM cols
        clinical_sub[missing_stages, "ajcc_pathologic_stage"] <- find_ajcc_stage(tnm_cols = tnm_sub)

        # omitt patinets you can't find survival data or stage data with
        missing_pat <- which(is.na(clinical_sub$OS_time) | is.na(clinical_sub$ajcc_pathologic_stage))
        # subsetting an empty list from the df will remove all entries, you must do a length check before removing
        if (length(missing_pat == 0)) clinical_sub <- clinical_sub[-missing_pat,]

        clinical_sub <- dplyr::select(clinical_sub, -c("ajcc_pathologic_t", "ajcc_pathologic_n", "ajcc_pathologic_m", "days_to_last_follow_up", "days_to_death", "days_to_diagnosis"))
    } else {
        clinical_sub <- dplyr::select(clinical_sub, -c("ajcc_pathologic_stage", "ajcc_pathologic_t", "ajcc_pathologic_n", "ajcc_pathologic_m", "days_to_last_follow_up", "days_to_death", "days_to_diagnosis"))
    }
    return(clinical_sub)
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

fetch_exp_data <- function(cancer_type, addBatch = TRUE, numericBatch = TRUE, scale = FALSE, primaryTumorOnly = FALSE, keepAllAliquot = FALSE, formatPatient = FALSE) {
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

    # Select the best aliquot
    if (!keepAllAliquot) {
        message(paste0("Expression matrix size went from :", dim(exp_matrix)[1], " * ", dim(exp_matrix)[2]))
        subset_bcr = filter_replicate_samples(rownames(exp_matrix))
        exp_matrix = exp_matrix[subset_bcr,]
        message(paste0("to :", dim(exp_matrix)[1], " * ", dim(exp_matrix)[2]))
    }
    
    # Convert all batch covaraites into numeric
    batches = c("TSS", "patient", "portion", "plate", "center")
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

    # Remove extra barcode entries for clinical data merge
    if(primaryTumorOnly & formatPatient) {
        rownames(exp_matrix) = format_tcga_patient(rownames(exp_matrix))
    } else if (!primaryTumorOnly & formatPatient) {
       warning("Formatting tumor barcode without first selecting for primary tumor will result in duplicated entries, not formatting patient barcode.")
    }

    return(exp_matrix)
}

mk_id_rownames <- function(df) {
    rownames(df) = df[,1]
    df = df[,-1]
    return(df)
}

convert_col_to_numeric <- function(clin_df, id = "donorId") {
    for (c in colnames(clin_df)) {
        if (!is.numeric(clin_df[,c]) && c != id) {
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
#' @param verbose print out which barcodes are kept and which ones are filtered
#' 
#' @return subset of the list of barcodes
#' 
#' @examples 
#' download.file("http://gdac.broadinstitute.org/runs/stddata__2014_01_15/samples_report/filteredSamples.2014_01_15__00_00_11.txt", "samples_filter.txt")
#'
#' test = read.table("samples_filter.txt", sep = "\t", header = T)
#' commas_bcr = unlist(strsplit(test_bcr[grepl(",", test_bcr)], ","))
#' test_bcr = c(test_bcr[!grepl(",", test_bcr)], commas_bcr, test$Chosen.Sample)
#' test_dna = test_bcr[grepl(".{19}[RHT]", test_bcr)]
#' test_rna = test_bcr[grepl(".{19}[DGWX]", test_bcr)]
#' filter_test_dna = filter_replicate_samples(test_dna, verbose = F)
#' filter_test_rna = filter_replicate_samples(test_rna, verbose = F)
#' if (all(filter_test_dna %in% test$Chosen.Sample) & all(filter_test_rna %in% test$Chosen.Sample)) { message("Filter passed.")
#' } else {message("Filter failed.")}
#' 
filter_replicate_samples <- function(bcr, verbose = TRUE) {
    if ( all(grepl(".{19}[RHT]", bcr)) ) { 
        type = "RNA"
    } else if (all (grepl(".{19}[DGWX]", bcr))) {
       type = "DNA"
    } else {
        stop("Mixing RNA and DNA samples not allowed.")
    }

    old_len <- length(bcr)

    bcr_df <- cbind(separate(as.data.frame(bcr),
                                    "bcr",
                                    c("project", "TSS", "patient", "sample,vial", "portion,analyte", "plate", "center"),  # skip var with NAs
                                    sep = "-"),
                                    bcr)

    bcr_df = separate(bcr_df, col = "sample,vial", into = c("sample", "vial"), sep = 2)
    bcr_df = separate(bcr_df, col = "portion,analyte", into = c("portion", "analyte"), sep = 2)

    # Selecting one aliquot based on GDC's Analyte Replicate Filter and Sort Replicate Filter rules
    bcr_df = bcr_df %>% arrange(analyte, desc(plate), desc(bcr)) %>% group_by(TSS, patient, sample) %>% dplyr::slice(1)
    
    new_len <- length(bcr_df$bcr)

    if (new_len == old_len) {
        print("No sample needs to be removed.")
        return(bcr_df$bcr)
    } else {
        out_tbl = data.frame()
        dup_samples = unique(substr(bcr[!(bcr %in% bcr_df$bcr)], 1, 15))
        for (dupbcr in dup_samples) {
            kept = as.character(bcr_df[grep(dupbcr, bcr_df$bcr),10])
            removed = bcr[grep(dupbcr, bcr)]
            removed = paste(removed[!(removed %in% kept)], collapse = ",")
            out_tbl = rbind(out_tbl, c(kept, removed))
        }
        colnames(out_tbl) = c("chosen", "removed")
        if (verbose) {
        print("The following changes are made: " )
        print(out_tbl)
        }
    }

    return(bcr_df$bcr)
}

#' Function for creating the treatment matrix (W) for each treatment variable. To pre-define treatment indicators before passing data into the `run.hte` function allows more flexibility.
#' 
#' @param txVector string vector indicating treatment variable names.
#' @param binaryVector bool vector indicating whether binary or continuous values were used. Default is set as binary for all treatment variables. Please remember to build this vector in concordance with the `is_binary` variable you have input in the `run.hte` function.
#' @param cutoffThreshDf numeric dataframe with first column indicating whether the group above this value is considered the treatment group or otherwise and the second column indicating the cutoff threshold. Default is taking entries with > 0.75 as the treatment group.
#' @param covarMat matrix containing per patient row entries and column entries of treamtent variable values.
#' 
#' 

create_tx_matrix <- function(txVector,
                            binaryVector = rep(TRUE, length(txVector)),
                            cutoffThreshDf = data.frame(
                                                dirct = rep(">", length(txVector)),
                                                thresh = rep(0.75, length(txVector))
                                                ),
                            covarMat = NULL) {
    W_matrix <- matrix(,nrow = dim(covarMat)[1], ncol = length(txVector))
    W_values <- NULL
    for (i in 1:length(txVector)) {
        W_values <- covarMat[,colnames(covarMat) == txVector[i]]
    	# print(W_values)
        if (binaryVector[i]){
           if (cutoffThreshDf[i, 1] == ">") {
               W_values <- as.numeric(W_values > quantile(W_values, cutoffThreshDf[i, 2]))
           } else if (cutoffThreshDf[i, 1] == "<") {
              W_values <- as.numeric(W_values < quantile(W_values, cutoffThreshDf[i, 2]))
           } else if (cutoffThreshDf[i, 1] == "0") {
              W_values <- as.numeric(W_values != 0)
           }
        }
	# print(W_values)
        W_matrix[, i] <- W_values
    }
    colnames(W_matrix) <- txVector
    return(W_matrix)
}

#' Correct drug names using a manually curated drug names conversion talbe provided by gatech.edu
#' You must first download the "DrugCorrection.csv" file to the ``./dat` directory from the gatech site via 
#' > wget https://gdisc.bme.gatech.edu/Data/DrugCorrection.csv --no-check-certificate
#' Unforutnately, setting `Sys.setenv(LIBCURL_BUILD="winssl")``, nor does setting `httr::set_config(config(ssl_verifypeer = FALSE))` work.
#' Update April 12, 2021, there are some capitalization and trailing space problems in the original table, use the corrected one "DrugCorrectionByAlex.csv"
#' 
#' @param drug_ls list of drugs extracted from TCGA datasets

correct_drug_names <- function(drug_ls) {
    drug_tbl = read.csv("./dat/DrugCorrectionByAlex.csv")
    drug_names <- c()
    for (i in 1:length(drug_ls)) {
        if (drug_ls[i] %in% drug_tbl$OldName) {
            drug_names[i] <- drug_tbl$Correction[drug_tbl$OldName == drug_ls[i]]
        } else {
            drug_names[i] <- drug_ls[i]
        }
    }
    return(drug_names)
}

#' Function to combine ajcc TNM system into stages according to the AJCC 7th edition guide, missing T/N/M indicators can be imputed
#' 
#' @example tnm_cols = dplyr::select(clinical_sub, all_of(c("ajcc_pathologic_t", "ajcc_pathologic_n", "ajcc_pathologic_m")))
#' @source \url{https://www.sciencedirect.com/science/article/pii/S1072751513002809}
#' @param tnm_cols subset of the clinical dataframe containing columns of pathological status of the TNM system
#' @param kSize k neighbors used for KNN impute
find_ajcc_stage <- function(tnm_cols = NULL, impute_missing = TRUE, kSize = 5) {

    message("Imputing missing stages.")
    if(impute_missing) {
        # impute missing T/N/M indicators by KNN
        tnm_impute <- sapply(tnm_cols, function(x) as.numeric(as.factor(x)))
        tnm_impute <- as.data.frame(knn.impute(tnm_impute, k = kSize))

        for (i in 1:ncol(tnm_cols)) {
            store_fac <- levels(factor(tnm_cols[,i]))
            store_num <- sort(unique(tnm_impute[,i]))
            store_convert <- as.data.frame(cbind(store_fac, store_num))
            
            # Convert the numeric factors from characters back to numeric
            store_convert$store_num = as.numeric(store_num)
            converted_factors <- left_join(data.frame(num_stage = tnm_impute[,i]), store_convert, by = c("num_stage" = "store_num")) 
            tnm_impute[,i] <- converted_factors$store_fac
        }
        tnm_cols <- tnm_impute
    }

    stage_tbl <- read.csv("./dat/brca_stage_tbl.csv")

    for (i in 1:ncol(tnm_cols)) {

        # remove minor categories
        with_subcat <- grep("^[[:upper:]]\\d[[:lower:]]$", tnm_cols[,i])
        tnm_cols[with_subcat,i] <- gsub("[[:lower:]]$", "", tnm_cols[with_subcat,i])

        # assume lowest possible stage
        with_unknown <- grep("^[[:upper:]]X$", tnm_cols[,i])
        tnm_cols[with_unknown,i] <- gsub("X$", "0", tnm_cols[with_unknown,i])

    }

    tnm_cols$key <- unite(tnm_cols, col = "key", sep = "|")
    stage_tbl$key <- unite(stage_tbl[,c(1:3)], col = "key", sep = "|")

    tnm_cols <- left_join(tnm_cols, stage_tbl, by = "key")

    # manually set for advanced stages
    tnm_cols[which(tnm_cols$ajcc_pathologic_n == "N3"), "stage"] <- "IIIC"
    tnm_cols[which(tnm_cols$ajcc_pathologic_m == "M1"), "stage"] <- "IV"

    return(tnm_cols$stage)
}

#' Function to fetch binary matrix indicating past pharmacological treatment taken by TCGA patients. Note that we don't assume this drug record is complete.
#' 
#' @source https://www.bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/clinical.html#Useful_information
#' @param cancerType TCGA project code
#' @param drugSummary whether to output number of people who has record of drug take, record of no drug take and no record at all
#' @return binary drug matrix of records avaiable on TCGA, not all patients are recorded on here though.
#' 
fetch_drug <- function(cancerType, drugSummary = TRUE) {

    message(paste0("Going through the TCGA BCR biotab records for: ", cancerType))
    query <- GDCquery(project = paste0("TCGA-", cancerType), 
                    data.category = "Clinical",
                    data.type = "Clinical Supplement", 
                    data.format = "BCR Biotab")
    GDCdownload(query)
    clinical.BCRtab.all <- GDCprepare(query)
    names(clinical.BCRtab.all)
    drug <- clinical.BCRtab.all[[paste0("clinical_drug_", tolower(cancerType))]]

    ## Clinical indexed data
    clinical <- GDCquery_clinic(project = paste0("TCGA-", cancerType), type = "clinical")

    # no drug record
    nr <- filter(clinical, treatments_pharmaceutical_treatment_or_therapy == "not reported")
    nr_drug <- drug[drug$bcr_patient_barcode %in% nr$bcr_patient_barcode,]

    # not taking any drug based on record
    no <- filter(clinical, treatments_pharmaceutical_treatment_or_therapy == "no")
    no_drug <- drug[drug$bcr_patient_barcode %in% no$bcr_patient_barcode,]

    # has at least one drug record
    yes <- filter(clinical, treatments_pharmaceutical_treatment_or_therapy == "yes")
    yes_drug <- drug[drug$bcr_patient_barcode %in% yes$bcr_patient_barcode,]

    if (drugSummary) {
        ## surveying the dataset
        message("In the official 'Clinical Indexed' data, among the patients with clinical entries, this is how much record we have: ")
        table(clinical$treatments_pharmaceutical_treatment_or_therapy)
        
        dim(nr)
            # [1] 103  74
        dim(nr_drug)
            # [1] 16 28
        message("Number of patients that was labeled as 'no record' in clinical indexed data but actually has drug record in the BCR biotab: ")
        message(length(unique(nr_drug$bcr_patient_barcode)), " out of ", length(unique(nr$bcr_patient_barcode)) )
            # 10 patient has drug record despite being labelled as "not reported"
        dim(no)
        # [1] 148  74
        dim(no_drug)
        # [1] 16 28
        message("Number of patients that was labeled as 'no drug taken' in clinical indexed data but actually has drug record in the BCR biotab: ")
        message(length(unique(no_drug$bcr_patient_barcode)), " out of ", length(unique(no$bcr_patient_barcode)) )
            # 5 patient has drug record despite being labelled as "no"
            # 143 patients are truly not taking drugs
        dim(yes)
        # [1] 846  74
        dim(yes_drug)
        # [1] 2378   28
        message("Number of patients that was labeled as 'yes' in clinical indexed data but actually has no drug record in the BCR biotab: ")
        message(length(unique(yes_drug$bcr_patient_barcode)), " out of ", length(unique(yes$bcr_patient_barcode)) )

        message("Corrected clinical indexed drug record: ")
        drug_summary <- c(rep("yes", length(unique(yes_drug$bcr_patient_barcode))), 
        rep("no", length(unique(no$bcr_patient_barcode)) - length(unique(no_drug$bcr_patient_barcode))), 
        rep("no record", length(unique(nr$bcr_patient_barcode)) - length(unique(nr_drug$bcr_patient_barcode))) )
        print(table(drug_summary))
    }
        # For each drug, we can choose the truely "no" patients as control, these are the universally non-treated patiets. If a patient with drug recrods is not treated with the drug in question but treated with other drugs in question, they will be the pseudo control group. We could try HTE with pseudo control or true control, depending on whichever works.

    # Fetch clinical indexed data from TCGAbiolinks
    clin_indexed <- fetch_indexed_clinical(cancerType)
    clin_indexed <- convert_col_to_numeric(clin_indexed, id = "bcr_patient_barcode")

    # Now the next step is to set up the drug binary tables, one hot encoding based on the above criteria

    drug_taken <- dplyr::select(drug, c("bcr_patient_barcode", "pharmaceutical_therapy_drug_name", "pharmaceutical_therapy_type"))
    drug_taken <- drug_taken[-c(1:2),]

    # Unify drug names with pre-downloaded conversion table
    drug_taken$corr_drug_names <- correct_drug_names(drug_taken$pharmaceutical_therapy_drug_name)
    drug_taken$corr_drug_names <- gsub(" ", "", drug_taken$corr_drug_names)
    drug_taken <- dplyr::select(drug_taken, -"pharmaceutical_therapy_drug_name") # get the old names out of the way

    ## Universal controls
    controls <- no$bcr_patient_barcode
    controls <- unique(controls[!(controls %in% drug_taken$bcr_patient_barcode)])

    ## True not reported cases (exclude)
    not_reported <- nr[!(nr$bcr_patient_barcode %in% drug_taken$bcr_patient_barcode),]
    not_reported <- not_reported$bcr_patient_barcode

    # co-prescription one hot encoding
    hot_drug <- dplyr::select(drug_taken, -c("pharmaceutical_therapy_type"))

    # re-enter multi-drug treated patients into multiple rows
    multi_tx <- hot_drug[grep("\\+", hot_drug$corr_drug_names),]
    if (nrow(multi_tx) != 0) { # only append multi drugs if there are patients taking multiple drugs, if not you will append an extra NA to the drug table
        for (i in 1:nrow(multi_tx)) {
            record <- unlist(strsplit(multi_tx$corr_drug_names[i], "\\+"))
            for (j in 1:length(record)) {
                hot_drug <- rbindlist(list(hot_drug, as.list(c(multi_tx$bcr_patient_barcode[i], record[j]))))
            }
        }
        hot_drug <- data.table(hot_drug[!grep("\\+", hot_drug$corr_drug_names)]) # removing rows by index not supported in data.tables
    }
    hot_drug <- dcast(setDT(hot_drug), bcr_patient_barcode ~ corr_drug_names, value.var = "corr_drug_names", function(x) as.numeric(length(x) > 0))

    # add true controls to same table
    control_row <- unique(no_drug$bcr_patient_barcode)

    for (ctrl in control_row) {
        hot_drug <- rbindlist(list(hot_drug, c(ctrl, as.list(rep(0, dim(hot_drug)[2] - 1)))))
    }
    return(list(hot_drug, clin_indexed, drug_taken, not_reported, controls))
}