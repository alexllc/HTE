
library(clusterProfiler)
# GO ontology enrichment

ensembl2ID <- AnnotationDbi::select(org.Hs.eg.db, keys=names(important),columns=c("ENTREZID"), keytype="ENSEMBL")

kk <- enrichKEGG(gene         = ensembl2ID$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)


## original sample
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]

kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)