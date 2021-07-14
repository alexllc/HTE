# Functions for retreiving detailed clinical data from all available sources (primarily from TCGAbiolinks)
# based on the TCGAbiolink documentation (https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/clinical.html), the XML files have more information, but the indexed data contains updated data with follow up information. And the BCRbiotabs are just the parsed tsv files from XML. So you will need to integrate both 

# tmp parameters:
cancer <- "KIRC"

# master function controlling input/output

#********************************** Cleaner functions ***********************************

unify_NAs <- function(df, na_words = c("Not Reported", "not reported", "NA", "<NA>", "NOS", "[Not Available]", "[Discrepancy]", "[Not Applicable]","[Not Available]","[Unknown]", "[Not Evaluated]")) {
    for(i in na_words) {
        df[df == i] <- NA
    }
    return(df)
}

dump_NA_cols <- function(df, ignore_rows = NULL, verbose = TRUE) {
    if (!is.null(ignore_rows))
        df <- df[-c(1:ignore_rows),]
    
    dump_list <- apply(df, 2, function(x) {all(is.na(x))})
    df <- df[,!dump_list]
    print("The following columns are discarded:")
    print(names(dump_list))
    return(df)
}

examine_df <- function(df) {

    for (col in colnames(df)) {
        tbl <- table(df[[col]])
        if(length(tbl) < 50) {
            message(paste(rep("=", 50), collapse = ""))
            message(col)
            print(tbl)
        } else {
            message(paste0("Variable skipped: ", col))
        }
    }
}

#********************************** downloading TCGA files ******************************

# Downloading BCR biotab
BCR_query <- GDCquery(project = paste0("TCGA-", cancer), 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab")
GDCdownload(BCR_query)
BCR_all <- GDCprepare(BCR_query)

# list all available tables
BCR_all_names <- names(BCR_all)

for (ds in BCR_all_names){
    tmp <- BCR_all[[ds]]
    tmp <- unify_NAs(tmp)
    tmp <- dump_NA_cols(tmp, ignore_rows = 2)
    assign(ds, tmp)
}


# Downloading indexed clincal files (a 75 col tbl)
indexed_clinical <- GDCquery_clinic(project = paste0("TCGA-", cancer), type = "clinical")

# clean and remove useless fields
indexed_clinical <- unify_NAs(indexed_clinical)
indexed_clinical <- dump_NA_cols(indexed_clinical, verbose = TRUE) # remove columns with only NAs entries

# summary columns
# smoker or non-smoker

# Retreive, re-calculate and impute (if necessary) tumor AJCC stage based on the TNM system

# general overview information
# stage, smoker, alcohol, 

# integrate other forms into the main indexed clinical table and provide means of queries if necessary
# clinical_drug_kirc: only add yes/no/no record to the table
# clinical_follow_up_v1.0_kirc: treatment_outcome_first_course, treatment_outcome_at_tcga_followup


# list to plot on graph
indexed_ls <- c(
    "gender",
    "race",
    "ajcc_pathologic_stage",
    "synchronous_malignancy",
    "prior_malignancy",
    "prior_treatment",
    "treatments_pharmaceutical_treatment_or_therapy",
    "treatments_radiation_treatment_or_therapy",
    "pack_years_smoked" # cigarette_per_year
)
sel_indexed <- dplyr::select(indexed_clinical, all_of(c("bcr_patient_barcode", indexed_ls)))

patient_ls <- c(
    "laterality",
    "race",
    "history_other_malignancy",
    "history_neoadjuvant_treatment",
    "lymph_nodes_examined_positive",
    "ajcc_pathologic_tumor_stage",
    "ldh_level",
    "tobacco_smoking_history_indicator",
    "radiation_treatment_adjuvant",
    "pharmaceutical_tx_adjuvant",
    "treatment_outcome_first_course",
    "new_tumor_event_dx_indicator"
)
sel_patient <- dplyr::select(clinical_patient_kirc, all_of(c("bcr_patient_barcode", patient_ls)))


follow_up_ls <- c(
    "radiation_treatment_adjuvant",
    "pharmaceutical_tx_adjuvant",
    "treatment_outcome_first_course",
    "treatment_outcome_at_tcga_followup",
    "new_tumor_event_dx_indicator",
    "new_tumor_event_surgery",
    "new_tumor_event_surgery_met",
    "new_tumor_event_radiation_tx",
    "new_tumor_event_pharmaceutical_tx"
)
sel_fu <- dplyr::select(clinical_follow_up_v1.0_kirc, all_of(c("bcr_patient_barcode", follow_up_ls)))

master_clinical <- full_join(sel_indexed, sel_patient, by = "bcr_patient_barcode")
master_clinical <- full_join(master_clinical, sel_fu, by = "bcr_patient_barcode")

# radiation details not required

# omf details not required yet

# nte details not required in plot yetq
indexed_clin_col <- colnames(master_clinical)[grep("\\.x$", colnames(master_clinical))]

for(col in indexed_clin_col) {
    var_name <- strsplit(col, "\\.")[[1]][1]
    indexed_tbl <- table(master_clinical[[col]])
    additional_tbl <- table(master_clinical[[paste0(var_name, ".y")]])
    names(indexed_tbl) <- toupper(names(indexed_tbl))
    names(additional_tbl) <- toupper(names(additional_tbl))
    if(!identical(indexed_tbl, additional_tbl)) {
        outdated <- ifelse(sum(additional_tbl) > sum(indexed_tbl), paste0(var_name, ".x"), paste0(var_name, ".y"))
    } else {
        outdated <- paste0(var_name, ".y")
    }
    master_clinical <- dplyr::select(master_clinical, -all_of(outdated))

}

col_suff <- grep("\\.x$|\\.y$", colnames(master_clinical))
tmp <- strsplit(colnames(master_clinical)[col_suff], "\\.")
tmp <- unlist(lapply(tmp, function(x) x[1][[1]]))
colnames(master_clinical)[col_suff] <- tmp

# Plot clinical data with tau values

# read tau values
tau <- read.csv("/home/alex/project/HTE/exp/2020-06-20_TCGA_pancancer_DEA_indicated_tx_dirct_.25_gp_as_tx_exp_HTE/res/KIRC/KIRC_tau_ENSG00000138180.csv")
tau <- tau[order(tau$tau.zval),]

# tau <- tau[!duplicated(tau$donorId),] # shouldn't need this if you didn't mess up with the data prep before forest
tau$direction <- ifelse(tau$tau.zval > 0, "Positive", "Negative")

pdf("tau_zscore_dot_ENSG00000138180.pdf")
ggdotchart(tau, x = "donorId", y = "tau.zval",
           color = "direction",                                # Color by groups
           palette = c("#E7B800", "#FC4E07"), # Custom color palette
           sorting = "ascending",                        # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           ggtheme = theme_pubr(),
           xlab = FALSE,
           ylab = "Tau Z score"                        # ggplot2 theme
           ) +
             theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                      axis.ticks.x=element_blank())
dev.off()