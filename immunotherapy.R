library(TCGAbiolinks)
library(DT) # for TCGAbiolinks
library(dplyr)
library(EDASeq)
library(tidyverse)

cancer_types = c(
    #'ACC','BLCA','BRCA',
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
usrwd = "/exeh_4/alex_lau/HTE/wd/expression_HTE"

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
    maf <- GDCquery_Maf(project, pipelines = "muse")

    # Below is the equivalence code
    # m_query <- GDCquery(project = "TCGA-HNSC",
    #               data.category = "Simple Nucleotide Variation",
    #               legacy = FALSE,
    #               workflow.type = "MuSE Variant Aggregation and Masking",
    #               data.type = "Masked Somatic Mutation"
    #             )

    # GDCdownload(m_query)
    # m_maf = GDCprepare(m_query)

    pmaf = dplyr::filter(maf, BIOTYPE == "protein_coding")
    pmaf$Hugo_Symbol[grep("ENSG00000085231", pmaf$Gene)] = "AK6"
    pmaf$kaks = ifelse(grepl("synonymous", pmaf$Consequence), "KS", "KA")
    gene_ln = as.data.frame(getGeneLengthAndGCContent(unique(pmaf$Gene), "hsa"))
    gene_ln$Gene = rownames(gene_ln)
    pmaf = left_join(pmaf, gene_ln, by = "Gene")
    pmaf = pmaf[!is.na(pmaf$length),] # some of the mut status has unknown gene

    # Create conversion table for later use
    hugo2ens = data.frame(Hugo_Symbol = as.character(pmaf$Hugo_Symbol), Gene = as.character(pmaf$Gene))
    ## ATTENTION: THE ORIGINAL MAF HAS A MISTAKE ON GENE NAME CONVERTION
    hugo2ens = hugo2ens[!duplicated(hugo2ens),]

    pmaf$SNNS = ifelse(pmaf$One_Consequence == "synonymous_variant" | is.na(pmaf$Amino_acids),"SN", "NS")
    # Using ==NA doesn't match with <NA>

    # Back ground mutation rate for EACH gene
    psn_by_gene = c()
    for (gene in unique(pmaf$Gene)) {
      n = gene_ln$length[gene_ln$Gene == gene]
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
    Pns_mat[is.na(Psn_mat)] = 0
    Pns_mat = cbind(wide$Tumor_Sample_Barcode, as.data.frame(Pns_mat))
    colnames(Psn_mat)[1] = "donorId"

    write.csv(Pns_mat, paste0(project, "_Pns.csv"), row.names=F)


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