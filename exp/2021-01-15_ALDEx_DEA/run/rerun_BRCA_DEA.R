# Tutorial see: https://bioc.ism.ac.jp/packages/3.2/bioc/vignettes/TCGAbiolinks/inst/doc/tcgaBiolinks.html#tcgaanalyze_dea-tcgaanalyze_leveltab-differential-expression-analysis-dea

library(TCGAbiolinks)

setwd("~/project/HTE/")

CancerProject <- "TCGA-BRCA"
DataDirectory <- paste0("./raw/",gsub("-","_",CancerProject))
FileNameData <- paste0(DataDirectory, "_","HTSeq_Counts",".rda")

query <- GDCquery(project = CancerProject,
                data.category = "Transcriptome Profiling",
                data.type = "Gene Expression Quantification", 
                workflow.type = "HTSeq - Counts")

samplesDown <- getResults(query,cols=c("cases"))

dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                typesample = "TP")

dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,
                                typesample = "NT")

queryDown <- GDCquery(project = CancerProject, 
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification", 
                    workflow.type = "HTSeq - Counts", 
                    barcode = c(dataSmTP, dataSmNT))
                    
GDCdownload(query = queryDown, directory = DataDirectory)

dataPrep <- GDCprepare(query = queryDown, directory =  DataDirectory)

# Filter out the samples with low correlations < 0.6 (outlier)
dataPrep <- TCGAanalyze_Preprocessing(object = dataPrep, 
                                    cor.cut = 0.6,
                                    datatype = "HTSeq - Counts")                      

dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                    geneInfo = geneInfoHT,
                                    method = "gcContent") 

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)   

# Diff.expr.analysis (DEA)
dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,dataSmNT],
                            mat2 = dataFilt[,dataSmTP],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")

# DEGs table with expression values in normal and tumor samples
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,"Tumor","Normal",
                                dataFilt[,dataSmTP],dataFilt[,dataSmNT])

write.csv(dataDEGsFiltLevel, file = "./dat/tables/BRCA_DEG_rerun.csv", row.names = TRUE)