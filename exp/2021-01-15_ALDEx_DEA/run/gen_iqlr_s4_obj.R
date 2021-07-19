setwd("~/project/HTE/")

library(ALDEx2)
library(TCGAbiolinks)
library(data.table)
library(gdata)
library(BiocParallel)
library(SummarizedExperiment)

# Fetch TCGA gene counts 
cancerList <- c("TCGA-KIRC", "TCGA-BLCA", "TCGA-COAD", "TCGA-GBM", "TCGA-HNSC", "TCGA-LUAD","TCGA-LUSC", "TCGA-PRAD", "TCGA-STAD", "TCGA-THCA", "TCGA-UCEC")

for (CancerProject in cancerList) {
    
    message(paste0("Creating IQLR object for cancer type:", CancerProject))

    if (!file.exists(paste0("./dat/iqlr_clr_obj/", CancerProject, "_iqlr_S4_obj.rds.gz"))) {
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

        if (!file.exists(paste0("./dat/HTSeq-count/TCGA_", strsplit(CancerProject, "-")[[1]][2], "_HTSeq_Counts.rda")) ) {

            # Donwload count matrix from GDC
            queryDown <- GDCquery(project = CancerProject, 
                                data.category = "Transcriptome Profiling",
                                data.type = "Gene Expression Quantification", 
                                workflow.type = "HTSeq - Counts", 
                                barcode = c(dataSmTP, dataSmNT))
                                
            GDCdownload(query = queryDown, directory = DataDirectory)

            data <- GDCprepare(query = queryDown, 
                            save = TRUE, 
                            save.filename = paste0("./raw/", CancerProject, "_count_matrix.RData"), # Ensembl server may go down 
                            directory =  DataDirectory)
            save(data, paste0("./raw/", CancerProject, "_count_matrix.RData"))

        } else {
            load(paste0("./dat/HTSeq-count/TCGA_", strsplit(CancerProject, "-")[[1]][2], "_HTSeq_Counts.rda"))
        }

        # convert HTSeq raw count into a data.frame
        cntmat <- data.frame(as.list(assays(data,withDimnames=TRUE)))

        # replace NAs with 0s
        cntmat[is.na(cntmat)] <- 0
        
        # Create condition list
        colnames(cntmat) <- gsub("HTSeq...Counts.", "", colnames(cntmat))
        conds <- gsub("\\.", "-", colnames(cntmat))
        names(conds) = ifelse(conds %in% dataSmTP, "TP", "NT")

        # Limit the number of cores available to ALDEx to avoid crashing
        default <- registered()
        register(MulticoreParam(workers = 40), default = TRUE)

        # unfiltered count matrix aldex
        message("Performing CLR.")
        iqlr <- aldex.clr(reads = cntmat, 
                        conds = names(conds), 
                        mc.samples = 128, 
                        denom="iqlr", 
                        verbose=TRUE, 
                        useMC=TRUE)

        saveRDS(iqlr, file = paste0("./dat/iqlr_clr_obj/", CancerProject, "_iqlr_S4_obj.rds.gz"))
        message("IQLR calculation done and saved.")
    } else {
        message("Reading IQLR object from previous run.")
        iqlr <- readRDS(paste0("./dat/iqlr_clr_obj/", CancerProject, "_iqlr_S4_obj.rds.gz"))
    }
    sampleIDs <- getSampleIDs(iqlr)
    expected_count <- matrix(1, nrow = numFeatures(iqlr), ncol = length(sampleIDs))

    for (ID in 1:length(sampleIDs)) {
        count <- rowMeans(iqlr@analysisData[[sampleIDs[ID]]])
        expected_count[, ID] <- count
    }
    rownames(expected_count) <- getFeatureNames(iqlr)
    colnames(expected_count) <- sampleIDs
    write.csv(expected_count, file = gzfile(paste0("./dat/iqlr_normalized_expr/", CancerProject, "_iqlr_expected_count.csv.gz")), row.names = TRUE)
    message("Point summary of IQLR object generated.")
}

