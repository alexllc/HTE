# Gene exp visualization

library(oncomix)

CancerProject <- paste0("TCGA-", cancer)
query <- GDCquery(project = CancerProject,
                data.category = "Transcriptome Profiling",
                data.type = "Gene Expression Quantification", 
                workflow.type = "HTSeq - Counts")

samplesDown <- getResults(query,cols=c("cases"))

dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                typesample = "TP")

dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,
                                typesample = "NT")

load(paste0("./HTSeqData/", project, "_exp.rda"))
prep = data
exp_matrix <- SummarizedExperiment::assay(prep, "HTSeq - FPKM-UQ")


mp = mixModelParams(exp_matrix[,dataSmNT], exp_matrix[,dataSmTP])
plot = plotGeneHist(mp, exp_matrix[,dataSmNT], exp_matrix[,dataSmTP], "ENSG00000003096")

library(ggplot2)
library(ggpubr)

a = ggplot(as.data.frame(tmp_tx), aes(x = tmp_tx))

a + geom_density(aes(y = ..count..), fill = "lightgray") +
  geom_vline(aes(xintercept = mean(tmp_tx)), 
             linetype = "dashed", size = 0.6,
             color = "#FC4E07")

ggsave("gph", plot = last_plot(), device = "pdf")