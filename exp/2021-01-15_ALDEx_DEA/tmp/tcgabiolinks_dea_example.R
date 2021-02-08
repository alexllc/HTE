# Script to perform DEA analysis by TCGAbiolinks

setwd("~/project/HTE/")

bin_ls = list.files("./bin")
for (bin in bin_ls){
    source(paste0("./bin/", bin))
}

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
dataSmTP_short <- dataSmTP[1:10]
dataSmNT_short <- dataSmNT[1:10]

queryDown <- GDCquery(project = CancerProject, 
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - Counts", 
                      barcode = c(dataSmTP_short, dataSmNT_short))
                    
GDCdownload(query = queryDown,
            directory = DataDirectory)

dataPrep <- GDCprepare(query = queryDown, 
                       save = TRUE, 
                       directory =  DataDirectory,
                       save.filename = FileNameData)

dataPrep <- TCGAanalyze_Preprocessing(object = dataPrep, 
                                      cor.cut = 0.6,
                                      datatype = "HTSeq - Counts")                      

dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfoHT,
                                      method = "gcContent") 

# boxplot(dataPrep, outline = FALSE)

# boxplot(dataNorm, outline = FALSE)

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)   

dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,dataSmTP_short],
                            mat2 = dataFilt[,dataSmNT_short],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")  

