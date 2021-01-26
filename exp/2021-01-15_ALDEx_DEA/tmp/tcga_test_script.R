setwd("~/project/HTE/")

library(ALDEx2)
library(TCGAbiolinks)
library(data.table)
library(gdata)

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

if (!file.exists("./dat/BRCA_post_GC_norm.csv.gz") ) {

    # Donwload count matrix from GDC
    queryDown <- GDCquery(project = CancerProject, 
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification", 
                        workflow.type = "HTSeq - Counts", 
                        barcode = c(dataSmTP, dataSmNT))
                        
    GDCdownload(query = queryDown, directory = DataDirectory)

    dataPrep <- GDCprepare(query = queryDown, 
                        save = FALSE, # Ensembl server may go down 
                        directory =  DataDirectory)

    # GC content normalization
    dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                        geneInfo = geneInfoHT,
                                        method = "gcContent")

    write.csv(dataNorm, file = gzfile("./dat/BRCA_post_GC_norm.csv.gz"))
} else {
    dataNorm <- read.csv("./dat/BRCA_post_GC_norm.csv.gz")
}

# Replace zeros with 1 according to propr's strategy
dataNorm[dataNorm == 0] <- 1

# Perform ALDEx2
rownames(dataNorm) <- dataNorm[, 1]
dataNorm <- dataNorm[, -1]
conds <- gsub("\\.", "-", colnames(dataNorm))
names(conds) = ifelse(conds %in% dataSmTP, "TP", "NT")

# unfiltered count matrix aldex
message("Performing CLR.")
iqlr <- aldex.clr(reads = dataNorm, 
                  conds = names(conds), 
                  mc.samples = 128, 
                  denom="iqlr", 
                  verbose=TRUE, 
                  useMC=TRUE)

saveRDS(iqlr, file = "./dat/BRCA_iqlr_S4_obj.rds.gz")

message("Done.")

# X <- aldex( reads = sel_dat,
#             conditions = names(conds),
#             mc.samples = 128,
#             test = "t",
#             effect = TRUE,
#             include.sample.summary = TRUE,
#             verbose = FALSE,
#             denom = "iqlr",
#             iterate = FALSE,
#             )


# TODO Comparison with previsouly calculated DE genes
# Agreement b/t ALDEx vs DEA vs CNA
