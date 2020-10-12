library(TCGAbiolinks)
library(dplyr)
  
    if (!file.exists(paste0("./Pns/TCGA-", project, "_Pns.csv"))) {
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
      Pns_mat = read.csv(paste0("./Pns/TCGA-", project, "_Pns.csv"))
    }