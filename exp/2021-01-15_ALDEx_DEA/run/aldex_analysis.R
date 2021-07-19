# Script to run the ALDEx2 package from an already existing CLR object

library(ALDEx2)
library(BiocParallel)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(TCGAbiolinks)

# Main function copied from the ALDEx2 package "aldex.r"
aldex_with_clr <- function(x, reads = NULL, conditions, mc.samples=128, test="t",
                  effect=TRUE, include.sample.summary=FALSE, verbose=FALSE, denom="all"){
  
  if(test == "iterative"){
    print("aldex.ttest: doing t-test")
    x.tt <- aldex.ttest(x, conditions, paired.test=FALSE)
    print("aldex.ttest: seeding a second t-test")
    nonDE.i <- which(rownames(reads) %in% rownames(x.tt[x.tt$wi.eBH > .05 | x.tt$we.eBH > .05, ]))
    if(length(nonDE.i) == 0) stop("no non-DE references found")
    x.tt <- aldex(reads, conditions, mc.samples=mc.samples, test="t",
                  effect=effect, include.sample.summary=include.sample.summary,
                  verbose=verbose, denom=nonDE.i)
  }else if(test == "t") {
    print("aldex.ttest: doing t-test")
    x.tt <- aldex.ttest(x, conditions, paired.test=FALSE)
  }else if(test == "glm"){
    print("aldex.glm: doing Kruskal Wallace and glm test")
    x.tt <- aldex.glm(x, conditions)
  }else{
    stop("argument 'test' not recognized")
  }
  
  if(effect == TRUE && test == "t"){
    print("aldex.effect: calculating effect sizes")
    x.effect <- aldex.effect(x, conditions,
                             include.sample.summary=include.sample.summary, verbose=verbose)
    z <- data.frame(x.effect, x.tt)
  }else{
    z <- data.frame(x.tt)
  }
  
  return(z)
}

setwd("~/proj/HTE/")

# Limit the number of cores available to ALDEx to avoid crashing
default <- registered()
register(MulticoreParam(workers = 50), default = TRUE)

# Load conditions
# Fetch TCGA gene counts
cancerList <- c("BLCA", "COAD", "GBM", "HNSC", "KIRC", "LUAD", "LUSC", "PRAD", "STAD", "THCA", "UCEC")

for (CancerProject in cancerList) {
  message(paste0("Processing cancer:", CancerProject))
  DataDirectory <- paste0("./raw/",gsub("-","_",CancerProject))
  FileNameData <- paste0(DataDirectory, "_","HTSeq_Counts",".rda")

  query <- GDCquery(project = paste0("TCGA-", CancerProject),
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")

  samplesDown <- getResults(query,cols=c("cases"))

  dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "TP")

  dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "NT")

  if (!file.exists(paste0("./raw/TCGA-", CancerProject, "_count_matrix.RData")) ) {

      # Donwload count matrix from GDC
      queryDown <- GDCquery(project = paste0("TCGA-", CancerProject),
                          data.category = "Transcriptome Profiling",
                          data.type = "Gene Expression Quantification", 
                          workflow.type = "HTSeq - Counts", 
                          barcode = c(dataSmTP, dataSmNT))
                          
      GDCdownload(query = queryDown, directory = DataDirectory)

      dataPrep <- GDCprepare(query = queryDown, 
                          save = TRUE, 
                          save.filename = paste0("./raw/", CancerProject, "_count_matrix.RData"), # Ensembl server may go down 
                          directory =  DataDirectory)

  } else {
      load(paste0("./raw/TCGA-", CancerProject, "_count_matrix.RData"))
  }

  # Create condition list
  cntmat <- data.frame(as.list(assays(data,withDimnames=TRUE)))
  colnames(cntmat) <- gsub("HTSeq...Counts.", "", colnames(cntmat))
  conds <- gsub("\\.", "-", colnames(cntmat))
  names(conds) = ifelse(conds %in% dataSmTP, "TP", "NT")

  message("Reading IQLR object from previous run.")
  x <- readRDS(paste0("./dat/TCGA-", CancerProject, "_iqlr_S4_obj.rds.gz"))

  message("Performing ALDEx2 with IQLR.")
  out <- aldex_with_clr(x = x, 
                        conditions = names(conds), 
                        include.sample.summary=TRUE, 
                        verbose=TRUE, 
                        denom="iqlr")

  saveRDS(out, file = paste0("./dat/", CancerProject, "_ALDEx2_result.rds.gz"))
  message("ALDEx2 done.")
}