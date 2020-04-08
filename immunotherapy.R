library(biomaRt)
library(TCGAbiolinks)
library(DT) # for TCGAbiolinks
library(dplyr)
library(EDASeq)
library(tidyverse)

cancer_types = c(
    # 'ACC','BLCA',
    'BRCA',
    'CESC','CHOL','COAD','DLBC','ESCA','GBM','HNSC','KICH', 'KIRC','KIRP','LGG','LIHC','LUAD','LUSC','MESO','OV','PAAD','PCPG','PRAD','READ', 'SKCM','STAD', 'TGCT','THCA','THYM','UCEC','UCS','UVM')

## Fixing drug names
# Doesn't work: download.file("https://gdisc.bme.gatech.edu/Data/DrugCorrection.csv", "DrugCorrection.csv")
#drug_tbl = read.csv("DrugCorrection.csv")

# correct_drug_names <- function(df) {
#   colnames(drug_tbl)[1] <- "drug_name"
#   df <- left_join(df, drug_tbl, by = "drug_name")
#   df$drug_name <- NULL
#   colnames(df)[colnames(df)=="Correction"] <- "drug_names"
#   return(df)
# }

# correct_drug_cat <- function(cat) {
#   targeted_rx <- read.csv(paste0(usrwd, "/cancerl/data/targeted_rx.csv"))
#   targeted_rx[,1] <- sapply(sapply(as.character(targeted_rx[,1]), function(x) strsplit(x, ' ')), function(x) x[[1]])
#   colnames(targeted_rx)[1] <- "drug_name"

# }
# usrwd = "/exeh_4/alex_lau/"
# setwd(paste0(usrwd, "/HTE/wd/immune_HTE"))
setwd("/exeh_4/alex_lau/HTE/wd/expression_HTE")

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


for (project in cancer_types) {
    print(paste0(paste0(rep('#', 20), collapse=""), "Downloading: ", project, paste0(rep('#', 20), collapse="")))
    #############################################################
    ### 1. Prepre TCGA clinical data
    #############################################################

    # Fetch survival data from GDC
    # query <- GDCquery(project = paste0("TCGA-", project),
    #             data.category = "Clinical",
    #             file.type = "xml") # Obtain list from TCGAbiolinks:::getProjectSummary(<enter_project>)

    # GDCdownload(query)
    # clinical_set <- c("drug", "patient")

    # for (c in clinical_set) {
    # assign(c, GDCprepare_clinic(query, clinical.info = c))
    # }

    # ### Select required information
    # s.drug = dplyr::select(drug, bcr_patient_barcode, days_to_drug_therapy_start,days_to_drug_therapy_end, drug_name, project)
    # ### Fix drug names
    # s.drug = correct_drug_names(s.drug)
    # s.drug = s.drug %>% mutate_if(is.factor, as.character) %>% mutate(duration = days_to_drug_therapy_end - days_to_drug_therapy_start) %>% select(-c(days_to_drug_therapy_end, days_to_drug_therapy_start))
    

    #############################################################
    ### 2. Prepare Expression data
    #############################################################

    # if (!file.exists(paste0("./HTSeqData/", project, "_exp.rda")) ) {

    # # Fetch expression data from GDC
    # g_query <- GDCquery(project = paste0("TCGA-", project),
    #                 data.category = "Transcriptome Profiling",
    #                 legacy = F,
    #                 data.type = "Gene Expression Quantification",
    #                 workflow.type = "HTSeq - FPKM-UQ",
    #                 sample.type = "Primary solid Tumor")

    # GDCdownload(g_query)
    # expdat <- GDCprepare(query = g_query,
    #                     save = TRUE,
    #                     save.filename = paste0("./HTSeqData/",project, "_exp.rda"))
    # }

    #############################################################
    ### 3. Prepare mutation data
    #############################################################
    # Seems like the simple one don't always work for all cancer types
    #maf <- GDCquery_Maf(project, pipelines = "muse")

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

    write.csv(Pns_mat, paste0("./Pns", project, "_Pns.csv"), row.names=F)


    #############################################################
    ### 4. Prepare immune proportion data
    #############################################################

    #############################################################
    ### 5. Prepare whole dataset,covariate matrix and tx vector
    #############################################################

} # next cancer


    #############################################################
    ### 6. Run HTE for all cancer types
    #############################################################