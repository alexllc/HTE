# Annotation, downloading, cleaning
library(biomaRt)
library(TCGAbiolinks)
library(DT) # for TCGAbiolinks
library(plyr)
library(dplyr)
library(EDASeq)
library(tidyverse)
library(readxl)

# Survival imputation
library(survival)

# For main function
library(doParallel)
library(grf)

library(MASS)
library(doMC)
library(data.table)
library(survminer)
library(doParallel)
library(methods)

# For causal forest
library(grf)
library(BART)
library(ranger)
library(randomForestSRC)
library(randomForest)


# GRF
setwd("/exeh_4/alex_lau/HTE/wd/")

source("./HTE/grf_parameters.R")
source("./HTE/HTE_main_functions.R")
source("./HTE/HTE_validation_functions.R")
source("./HTE/survival_imputation.R")

cancer_types = c(
    # 'ACC','BLCA','BRCA', 'CESC','CHOL','COAD','DLBC','ESCA',
    'GBM',
    'HNSC','KICH', 'KIRC','KIRP','LGG','LIHC','LUAD','LUSC','MESO','OV','PAAD','PCPG','PRAD','READ', 'SKCM','STAD', 'TGCT','THCA','THYM','UCEC','UCS','UVM')

## Fixing drug names
# Doesn't work: download.file("https://gdisc.bme.gatech.edu/Data/DrugCorrection.csv", "DrugCorrection.csv")
drug_tbl = read.csv("./immunotx/DrugCorrection1.csv")

correct_drug_names <- function(df) {
  colnames(drug_tbl)[1] <- "drug_name"
  df <- left_join(df, drug_tbl, by = "drug_name")
  df$drug_name <- NULL
  colnames(df)[colnames(df)=="Correction"] <- "drug_names"
  return(df)
}

targeted_rx <- read.csv("./immunotx/targeted_rx.csv")
correct_drug_cat <- function(cat) {
  targeted_rx[,1] <- sapply(sapply(as.character(targeted_rx[,1]), function(x) strsplit(x, ' ')), function(x) x[[1]])
  colnames(targeted_rx)[1] <- "drug_name"

}

targetedrx = read.csv("./immunotx/targeted_rx.csv")
targetedrx = strsplit(as.character(targetedrx$Agent), ' ')
targetedrx = unlist(lapply(targetedrx, function(x) x[[1]]))
tcgatargrx = targetedrx[targetedrx %in% drug_tbl$Correction]

# usrwd = "/exeh_4/alex_lau/"
# setwd(paste0(usrwd, "/HTE/wd/immune_HTE"))

# Faster way of getting gene length?

# if(!file.exists) {
#   ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
#   proteinESG = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filter = c("transcript_biotype"), values = c("protein_coding"), mart = ensembl)
#   prot_sub = proteinESG[-grep("MT-*", proteinESG$hgnc_symbol),]
#   full_protls = unique(prot_sub$ensembl_gene_id)
#   geneLN1 = as.data.frame(getGeneLengthAndGCContent(unique(full_protls[1:500]), "hsa"))
#   geneLN2 = as.data.frame(getGeneLengthAndGCContent(unique(full_protls[1:2500]), "hsa"))
#   write.csv(geneLN, "geneLN.csv", row.names = F)
# }

# That didn't work for some reasons
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# library(org.Hs.eg.db)
library(Homo.sapiens)
keytypes(org.Hs.eg.db)
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
txlen = transcriptLengths(txdb, with.utr5_len=F, with.utr3_len=F)
txnames = AnnotationDbi::select(Homo.sapiens, txlen$tx_name, c("ENSEMBL", "SYMBOL"), "TXNAME")
txlen = left_join(txlen, txnames, by = c("tx_name" = "TXNAME"))
# txnames = AnnotationDbi::select(org.Hs.eg.db, unique(na.omit(txlen$gene_id)), c("ENSEMBL", "SYMBOL"), "ENTREZID")
txlen = dplyr::select(txlen, c(ENSEMBL, SYMBOL, tx_len))
txlen = txlen[complete.cases(txlen),]
txlen = txlen %>% group_by(SYMBOL) %>% mutate(avglen = mean(tx_len)) %>% ungroup() %>% group_by(ENSEMBL) %>% dplyr::slice((1))
txlen$tx_len = NULL
txlen = txlen[!duplicated(txlen),]
# Those that are protein coding will have their own ENSEMBL id, so selecting for entries with ENSEMBL != NA automatically selects for protein coding genes only


for (project in cancer_types) {
    print(paste0(paste0(rep('#', 20), collapse=""), "Downloading: ", project, paste0(rep('#', 20), collapse="")))
    #############################################################
    ### 1. Prepre TCGA clinical data
    #############################################################

    if (!file.exists(basename("TCGA_CDR_clean.csv"))) {

      download.file("https://ars.els-cdn.com/content/image/1-s2.0-S0092867418302290-mmc1.xlsx", "TCGA-CDR-SupplementalTableS1.xlsx")

      cdr = read_excel("TCGA-CDR-SupplementalTableS1.xlsx")
      clinical_dat = dplyr::select(cdr, c(bcr_patient_barcode, type, age_at_initial_pathologic_diagnosis,  gender, ajcc_pathologic_tumor_stage, OS, OS.time))

      #patient <- patient[patient$OS.time!=0,]
      clinical_dat[clinical_dat == "#N/A"] <- NA
      clinical_dat <- subset(clinical_dat, !is.na(OS.time) & !is.na(OS) & !is.na(age_at_initial_pathologic_diagnosis))
      colnames(clinical_dat)[colnames(clinical_dat)=="bcr_patient_barcode"] = "donorId"
      # This renaming step is critical as the HTE main function will rely on the column named outcome to indicate Y
        #patient$gender = as.numeric(patient$gender == 'FEMALE')
        write.csv(clinical_dat, "TCGA_CDR_clean.csv", row.names = F)
    } else {
        clinical_dat = read.csv("TCGA_CDR_clean.csv")
    }

    ss_patient <- subset(clinical_dat, type %in% project & OS.time > 0)
    surv.times <- as.numeric(ss_patient$OS.time)
    cens <- as.numeric(ss_patient$OS)

    # get imputed log survival times
    max.censored <- max(surv.times[cens == 0])
    cens[surv.times == max.censored] <- 1
    outcome = impute.survival(surv.times, cens)

    # attach imputed.log.times to original dataset
    ss_patient <- cbind(ss_patient, outcome)

    for (c in colnames(ss_patient)) {

        if (!is.numeric(ss_patient[,c]) && c != "donorId") {
            which.one <- which( levels(ss_patient[,c]) == "")
            levels(ss_patient[,c])[which.one] <- NA
            ss_patient[,c] = sapply(sapply(ss_patient[,c], as.factor), as.numeric) 
            print(paste0(c, " is altered")) 
        }
    }
    ss_patient = dplyr::select(ss_patient, -c(type, OS, OS.time))
    print("Processed patient dataframe: ")
    head(ss_patient)

    ### We will attach survival last

    # Fetch clinical data from GDC
    query <- GDCquery(project = paste0("TCGA-", project),
                data.category = "Clinical",
                file.type = "xml") # Obtain list from TCGAbiolinks:::getProjectSummary(<enter_project>)

    GDCdownload(query)
    clinical_set <- c("drug", "patient")

    for (c in clinical_set) {
    assign(c, GDCprepare_clinic(query, clinical.info = c))
    }

    ### Select required information
    s.drug = dplyr::select(drug, bcr_patient_barcode, drug_name, project)
    # days_to_drug_therapy_start,days_to_drug_therapy_end, drug_name, project)
    # ### Fix drug names
    # s.drug = correct_drug_names(s.drug)
    # s.drug = s.drug %>% mutate_if(is.factor, as.character) %>% mutate(duration = days_to_drug_therapy_end - days_to_drug_therapy_start) %>% dplyr::select(-c(days_to_drug_therapy_end, days_to_drug_therapy_start))
    s.drug  = s.drug[!is.na(s.drug$drug_name),]
    colnames(s.drug)[1] = "donorId"
    s.drug = s.drug %>% group_by(donorId) %>% mutate(chemo_only = !(any(drug_name %in% tcgatargrx)))
    for (d in tcgatargrx) {
      s.drug[d] = s.drug$drug_name == d
    }
    s.drug = s.drug %>% group_by(donorId) %>% mutate_if(is.logical, any) %>% mutate_if(is.logical, as.numeric)
    s.drug$drug_name = NULL
    s.drug = s.drug[!duplicated(s.drug),]

    clinical = left_join(s.drug, ss_patient, by = "donorId")

    #############################################################
    ### 2. Prepare mutation data
    #############################################################
    # Seems like the simple one don't always work for all cancer types
    #maf <- GDCquery_Maf(project, pipelines = "muse")

    if (!file.exists(paste0("./Pns/TCGA-", project, "_Pns.csv"))){
    # Below is the equivalence code
    m_query <- GDCquery(project = paste0("TCGA-", project),
                  data.category = "Simple Nucleotide Variation",
                  legacy = FALSE,
                  workflow.type = "MuSE Variant Aggregation and Masking",
                  data.type = "Masked Somatic Mutation"
                )

    GDCdownload(m_query)
    maf = GDCprepare(m_query)

    pmaf = dplyr::filter(maf, BIOTYPE == "protein_coding")
    # pmaf$Hugo_Symbol[grep("ENSG00000085231", pmaf$Gene)] = "AK6"
    pmaf$kaks = ifelse(grepl("synonymous", pmaf$Consequence), "KS", "KA")
    pmaf = left_join(pmaf, txlen, by = c("Gene" = "ENSEMBL"))
    pmaf = pmaf[!is.na(pmaf$avglen),] # some of the mut status has unknown gene

    # Create conversion table for later use
    hugo2ens = data.frame(Hugo_Symbol = as.character(pmaf$SYMBOL.y), Gene = as.character(pmaf$Gene))
    ## ATTENTION: THE ORIGINAL MAF HAS A MISTAKE ON GENE NAME CONVERTION
    hugo2ens = hugo2ens[!duplicated(hugo2ens),]

    pmaf$SNNS = ifelse(pmaf$One_Consequence == "synonymous_variant" | is.na(pmaf$Amino_acids),"SN", "NS")
    # Using ==NA doesn't match with <NA>

    # Back ground mutation rate for EACH gene
    psn_by_gene = c()
    for (gene in unique(pmaf$Gene)) {
      n = floor(txlen$avglen[txlen$ENSEMBL == gene])
      x = sum(pmaf$Gene == gene & pmaf$SNNS == "SN")
      psn = binom.test(x,n)
      psn = as.numeric(psn$estimate)
      psn = c(gene, psn)
      psn_by_gene = rbind(psn_by_gene, psn)
    }
    psn_by_gene = as.data.frame(psn_by_gene)
    colnames(psn_by_gene) = c("Gene", "Psn")
    #rownames(psn_by_gene) = NULL
    psn_by_gene = left_join(psn_by_gene, hugo2ens, by = "Gene")
    psn_by_gene = psn_by_gene[!duplicated(psn_by_gene),]
    psn_by_gene = psn_by_gene[order(psn_by_gene$Hugo_Symbol),]

    # Calcualte Pns for each patient each gene
    pmaf = pmaf %>% group_by(Tumor_Sample_Barcode) %>% mutate(R = sum(SNNS == "NS") / sum(SNNS == "SN"))
    # Some of the patients only have NS therefore R becomes infinity, we will give them the second largest value
    pmaf$R[pmaf$R == Inf] = sort(unique(pmaf$R), decreasing=T)[2]

    # subset and convert to wide format
    # We will use Hugo symbol for mutation while ens is reserved for expression
    spmaf = dplyr::select(pmaf, c(Hugo_Symbol, Tumor_Sample_Barcode, R))
    spmaf = spmaf[!duplicated(spmaf),]
    wide = spmaf %>% spread(Hugo_Symbol, R)

    # Get Pns by multiplying R with the Psn of EACH gene
    Pns = as.numeric(as.character(psn_by_gene$Psn))
    tmp = as.matrix(wide[,2:ncol(wide)])
    class(tmp) = "numeric"
    Pns_mat = t(t(tmp) * Pns)
    Pns_mat[is.na(Pns_mat)] = 0
    Pns_mat = cbind(wide$Tumor_Sample_Barcode, as.data.frame(Pns_mat))
    colnames(Pns_mat)[1] = "donorId"
    tmp = strsplit(as.character(Pns_mat$donorId), "-")
    tmp = unlist(lapply(tmp, function(x) paste(x[[1]], x[[2]], x[[3]], sep = "-")))
    Pns_mat$donorId <- unlist(tmp)
    } else{
      Pns_mat = paste0("./Pns/TCGA-", project, "_Pns.csv"))
    }
    

    #write.csv(Pns_mat, paste0("./Pns/", project, "_Pns.csv"), row.names=F)


    #############################################################
    ### 3. Prepare immune proportion data
    #############################################################
    immuneprop = read.csv(paste0("./immunotx/proportion/", project, "_immune_cells.csv"))


    #############################################################
    ### 4. Prepare Expression data
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
    } else {
        load(paste0("./HTSeqData/", project, "_exp.rda"))
        prep = data
    }

    exp_matrix <- SummarizedExperiment::assay(prep, "HTSeq - FPKM-UQ")

    # Only select primary tumor samples
    # Somehow this line removed so many patients
    exp_matrix <- exp_matrix[,grep("-01.-", colnames(exp_matrix))]

    #Scale FPKM value by patients
    exp_matrix <- exp_matrix/apply(exp_matrix,2,max)
    exp_matrix <- t(exp_matrix)

    exp_matrix = dplyr::select(as.data.frame(exp_matrix), colnames(exp_matrix) [colnames(exp_matrix) %in% txlen$ENSEMBL])

    #Correct for batch effects
    exp_matrix <- cbind(separate(as.data.frame(rownames(exp_matrix)), "rownames(exp_matrix)",
                            c(NA, "TSS", "patient", NA, "portion", "plate", "center"),  # skip var with NAs
                            sep = "-"), exp_matrix)
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


    #############################################################
    ### 5. Prepare whole dataset,covariate matrix and tx vector
    #############################################################
    wholedat = left_join(clinical, Pns_mat, by = "donorId") %>% left_join(. , exp_matrix, by = "donorId")
    wholedat$project = as.numeric(as.factor(wholedat$project))
    assign(paste0("wholedat", project), wholedat)
    write.csv(get(paste0("wholedat", project)), paste0("./immunotx/wholedats/wholedat_", project), row.names = F)
    print(paste0(project, " whole dataset saved."))

} # next cancer
    
    # pancan_whole = as.data.frame(rbind.fill(wholedatLUAD, wholedatLUSC)) # for combining all patients
    # # quickest way known to R for making NA into 0
    # for (j in seq_len(ncol(pancan_whole))) {
    #   set(pancan_whole,which(is.na(pancan_whole[[j]])),j,0)
    # }
    # pancan_whole$project = as.numeric(as.factor(pancan_whole$project))
    # #############################################################
    # ### 6. Run HTE for all cancer types
    # #############################################################
    # txdrug = "Gefitinib"
    # covar_mat= as.data.frame(dplyr::select(pancan_whole, -c("donorId", "outcome")))

    # # Run HTE with the whole dataset and covar matrix
    # obsNumber <- dim(covar_mat)[1]
    # trainId <- sample(1: obsNumber, floor(obsNumber/2), replace = FALSE)
    # registerDoParallel(10)
    # output_file = "./test_output/"
    # result <- run.hte(covar_mat, txdrug, pancan_whole, project, covar_type = "drug", trainId, seed = 111, is.binary = F, is_save = T, save_split = T, is.tuned = F, thres = 0.75, n_core = 8, output_directory = output_file)
    # write.csv(result[[1]], paste0(output_file, project, '_expression_correlation_test_result.csv'), quote = F, row.names = F)
    # write.csv(result[[2]], paste0(output_file, project, '_expression_calibration_result.csv'), quote = F, row.names = F)
    # write.csv(result[[3]], paste0(output_file, project, '_expression_median_t_test_result.csv'), quote = F, row.names = F)
    # write.csv(result[[4]], paste0(output_file, project, '_expression_permutate_testing_result.csv'), quote = F, row.names = F)