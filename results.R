library(ReactomePA)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(Homo.sapiens)
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene

cancer_list = c(
                'BLCA',
                'COAD',
                'BRCA',
                #'LGG',
                'GBM',
                'STAD',
                'HNSC',
                'KIRC',
                'LUAD',
                'LUSC',
                #'OV', # no NT
                'PRAD',
                #'SKCM', # no NT
                'THCA',
                'UCEC',
                'ESCA')

total = 0
for (c in cancer_list) {
    tmp = read.csv(paste0("./", c, "/", c, "_expression_permutate_testing_result.csv"))
    tmp = tmp[tmp[,2] < 0.05 & tmp[,3] < 0.05,]
    assign(paste0(c, ".perm"), tmp)
    total = total + dim(tmp)[1]
    print(paste0(c, " has ", dim(tmp)[1], " HTE genes."))
}

txnames = AnnotationDbi::select(Homo.sapiens, as.character(varimp$variable), c("ENTREZID", "SYMBOL"), "ENSEMBL")
x <- enrichPathway(gene=txnames$ENTREZID,pvalueCutoff=0.05, readable=T)
head(as.data.frame(x))

stats = NULL
for (gene in BRCA.perm$mutName){
    tmp = read.csv(paste0("./", c, "/", c, "_tau_", gene, ".csv"))
    dat = c(max(tmp$tau.val), min(tmp$tau.val), abs(exp(max(tmp$tau.val)) - exp(min(tmp$tau.val))))
    stats = rbind(stats, dat)
    
}
stats = as.data.frame(stats)
stats = stats[order(stats$V3),]

