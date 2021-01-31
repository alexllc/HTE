setwd("~/proj/HTE/")

library(ALDEx2)
library(TCGAbiolinks)
library(data.table)
library(gdata)
library(BiocParallel)
library(SummarizedExperiment)

# Fetch TCGA gene counts
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

if (!file.exists("./raw/count_matrix.RData") ) {

    # Donwload count matrix from GDC
    queryDown <- GDCquery(project = CancerProject, 
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification", 
                        workflow.type = "HTSeq - Counts", 
                        barcode = c(dataSmTP, dataSmNT))
                        
    GDCdownload(query = queryDown, directory = DataDirectory)

    dataPrep <- GDCprepare(query = queryDown, 
                        save = TRUE, 
                        save.filename = "./raw/count_matrix.RData", # Ensembl server may go down 
                        directory =  DataDirectory)

    } else {
    load("./raw/count_matrix.RData")
}

# Create condition list
cntmat <- data.frame(as.list(assays(data,withDimnames=TRUE)))
colnames(cntmat) <- gsub("HTSeq...Counts.", "", colnames(cntmat))
conds <- gsub("\\.", "-", colnames(cntmat))
names(conds) = ifelse(conds %in% dataSmTP, "TP", "NT")

# Limit the number of cores available to ALDEx to avoid crashing
default <- registered()
register(MulticoreParam(workers = 28), default = TRUE)

# unfiltered count matrix aldex
message("Performing CLR.")
iqlr <- aldex.clr(reads = cntmat, 
                  conds = names(conds), 
                  mc.samples = 128, 
                  denom="iqlr", 
                  verbose=TRUE, 
                  useMC=TRUE)

saveRDS(iqlr, file = "./dat/BRCA_iqlr_S4_obj.rds.gz")
message("IQLR calculation done and saved.")

# Summarize MC instances into one expected value
for ()


