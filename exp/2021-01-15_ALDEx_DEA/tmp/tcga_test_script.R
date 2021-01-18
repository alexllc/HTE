setwd("~/project/HTE/")

library(ALDEx2)
library(TCGAbiolinks)

# Fetch TCGA gene counts
CancerProject <- "TCGA-BRCA"
DataDirectory <- "./exp/2021-01-15_ALDEx_DEA/raw"

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

dataPrep <- GDCprepare(query = queryDown, 
                    save = TRUE, 
                    directory =  DataDirectory)

# Filter out zeros and near zero values
dataFilt <- TCGAanalyze_Filtering(tabDF = dataPrep,
                                method = "quantile", 
                                qnt.cut =  0.25)
# Perform ALDEx2
conds <- colnames(dataPrep)
names(conds) = ifelse(conds %in% dataSmTP, "TP", "NT")
conds_filt <- colnames(dataFilt)
names(conds_filt) = ifelse(conds_filt %in% dataSmTP, "TP", "NT")

# TODO unfiltered count matrix aldex

iqlr <- aldex.clr(reads = dataPrep, conds = conds, mc.samples = 128, denom="iqlr", verbose=FALSE, useMC=TRUE)


X <- aldex( reads = dataPrep,
            conditions = conds,
            mc.samples = 128,
            test = "t",
            effect = TRUE,
            include.sample.summary = TRUE,
            verbose = FALSE,
            denom = "iqlr",
            iterate = FALSE,
            )


# TODO Comparison with previsouly calculated DE genes
