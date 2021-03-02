## Run from base directory

setwd("~/project/HTE/")
bin_ls = list.files("./bin")
for (bin in bin_ls){
    source(paste0("./bin/", bin))
}

library(org.Hs.eg.db)
library(dplyr)
library("RColorBrewer")

get_gene <- function(file_name, pos = NULL) {
    gene_names <- strsplit(file_name, "_|\\.")
    gene_names <- unlist(lapply(gene_names, function(x) {x[pos]}))
    return(gene_names)
}


# heatmap parameters
expri <- "/home/alex/project/HTE/exp/2020-06-20_TCGA_pancancer_DEA_indicated_tx_dirct_.25_gp_as_tx_exp_HTE"
cancer <- "BRCA"
include_rx <- TRUE

perm <- read.csv(paste0(expri, "/res/", cancer, "/", cancer, "_expression_permutate_testing_result.csv"))
# Load perm results and filter tau by it
sig_perm <- filter(perm, var.pval < 0.05 | fixed.YW.risk.pval < 0.05)

# Load tau files
tau_files <- list.files(path = paste0(expri, "/res/", cancer, "/"), pattern = "_tau_")
tau_tbl <- data.frame()
for (i in 1:length(tau_files)) {
    gene_tau <- tau_files[i]
    tx_gene <- get_gene(gene_tau, pos = 3)
    tau_values <- read.csv(paste0(expri, "/res/", cancer, "/", gene_tau))
    tau_vec <- tau_values$tau.val
    if (i == 1) names(tau_vec) <- tau_values$donorId
    tau_tbl <- rbind(tau_tbl, as.list(tau_vec))
}

# Convert ENSEMBL gene names to Hugo Symbols
ensembl_gene <- get_gene(tau_files, pos = 3)
db_translate <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys=ensembl_gene,columns=c("SYMBOL"), keytype="ENSEMBL"))

# Sometimes one SYMBOL maps to multiple transcripts
symbol_unique <- NULL
for (gene in unique(db_translate$ENSEMBL)) {
    symbol_ls <- db_translate[db_translate$ENSEMBL == gene,]
    symbol_ls <- paste(symbol_ls$SYMBOL, collapse = "|")
    symbol_unique <- c(symbol_unique, symbol_ls)
}
rownames(tau_tbl) <- ensembl_gene
tau_tbl <- tau_tbl[rownames(tau_tbl) %in% sig_perm$gene,]
colnames(tau_tbl) <- gsub("\\.", "-", colnames(tau_tbl))

# include important clinical infomration
cdr <- read_excel("./dat/TCGA-CDR-SupplementalTableS1.xlsx")
clinical <- dplyr::select(cdr, all_of(c("bcr_patient_barcode", "age_at_initial_pathologic_diagnosis", "gender", "ajcc_pathologic_tumor_stage", "tumor_status")))
clinical <- clinical[clinical$bcr_patient_barcode %in% colnames(tau_tbl),]

# Load drug taken from TCGA BCR Biotab files
query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab")
GDCdownload(query)
clinical.BCRtab.all <- GDCprepare(query)
names(clinical.BCRtab.all)
drug <- clinical.BCRtab.all$clinical_drug_brca
drug_taken <- dplyr::select(drug, c("bcr_patient_barcode", "pharmaceutical_therapy_drug_name", "pharmaceutical_therapy_type"))
drug_taken <- drug_taken[-c(1:2),]

# Unify drug names with pre-downloaded conversion table
drug_taken$corr_drug_names <- correct_drug_names(drug_taken$pharmaceutical_therapy_drug_name)
drug_taken <- dplyr::select(drug_taken, -"pharmaceutical_therapy_drug_name") # get the old names out of the way

# We only care about whether a patient has *taken* the drug, not for how long or how many regimes, this information is important but it adds extra complexity in our visualization
drug_taken <- dplyr::select(drug_taken, -"pharmaceutical_therapy_type")
# long to wide drug information
# Unique drug entries for each patients
wide_drug <- drug_taken %>% unique() %>% mutate(n = 1) %>% spread(corr_drug_names, n, fill = 0)

#************************ build heat map ************************

plot_clinical <- data.frame(stage = clinical$ajcc_pathologic_tumor_stage, 
                        gender = clinical$gender, 
                        tumor_status = clinical$tumor_status,
                        age = as.character(cut(clinical$age_at_initial_pathologic_diagnosis, c(seq(20, 90, by = 10), include.lowest = TRUE)))
                        )

# replace NAs with a character so pheatmap won't reject our plot
plot_clinical[is.na(plot_clinical)] <- "Not available"

rownames(plot_clinical) <- clinical$bcr_patient_barcode

stage <- colorRampPalette(brewer.pal(length(unique(plot_clinical$stage)) - 2, "PuRd"))(length(unique(plot_clinical$stage)) - 2)
names(stage) <- sort(unique(plot_clinical$stage))[3:length(unique(plot_clinical$stage))]
stage <- c(stage, c("[Not Available]" = "azure4", "[Discrepancy]" = "lemonchiffon4"))

age <- brewer.pal(length(unique(plot_clinical$age)), "BuGn")
names(age) <- sort(unique(plot_clinical$age))

quali_col <- list(stage = stage,
                  gender = c("FEMALE" = "gold1", "MALE" = "blue4"),
                  tumor_status = c("WITH TUMOR" = "#DE2D26", "TUMOR FREE" = "#FEE0D2", "Not available" = "honeydew4"),
                  age = age
                )

# include drug tx categories
if (include_rx) {
    plot_clinical_rx <- plot_clinical
    plot_clinical_rx$bcr_patient_barcode <- rownames(plot_clinical_rx)
    plot_clin_rx <- left_join(plot_clinical_rx, wide_drug, by = "bcr_patient_barcode")
    rownames(plot_clin_rx) = plot_clin_rx$bcr_patient_barcode
    plot_clin_rx$bcr_patient_barcode = NULL
    plot_clin_rx[is.na(plot_clin_rx)] = "missing" # do not indicate missing drug information as "0", not taken
    plot_clinical <- plot_clin_rx

    drug_col <- list()
    for (i in 5:dim(plot_clin_rx)[2]) {
        x <- c("1" = "brown1", "0" = "gray88", "missing" = "white")
        assign(colnames(plot_clin_rx)[i], x)
        drug_col[[i-4]] <- get(colnames(plot_clin_rx)[i])
    }
    names(drug_col) <- colnames(plot_clin_rx)[5:dim(plot_clin_rx)[2]]
    quali_col <- append(quali_col, drug_col)
}

qualitative_var <- list(stage <- clinical$ajcc_pathologic_tumor_stage, 
                        gender <- clinical$gender, 
                        tumor_status <- clinical$tumor_status,
                        age <- cut(clinical$age_at_initial_pathologic_diagnosis, c(seq(20, 90, by = 10), include.lowest = TRUE)),
                        col = quali_col
                        )

pdf(file = "HTE_heatmap.pdf", width = 30, height = 20)
pheatmap(as.matrix(tau_tbl),
        annotation_col = plot_clinical,
        annotation_colors = quali_col,
        annotation_legend = F,
        cluster_cols=T, cluster_rows=T,
        scale = "row",
        show_colnames = F
        )
dev.off()