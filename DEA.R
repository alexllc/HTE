library(TCGAbiolinks)
library(SGSeq)
library(biomaRt)

setwd("/home/alex/project/HTE/wd/")

cancer_list = c(
                # 'BLCA',
                'COAD',
                'BRCA',
                # 'LGG', # no NT
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

#ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

for (cancer in cancer_list) {
    
    print(paste0("Begin processing: ", cancer))
    CancerProject <- paste0("TCGA-", cancer)
    DataDirectory <- paste0("./expression_HTE/DEA/",gsub("-","_",CancerProject))
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

    dataPrep <- GDCprepare(query = queryDown, 
                        save = TRUE, 
                        directory =  DataDirectory,
                        save.filename = FileNameData)
    # Filter out the samples with low correlations < 0.6 (outlier)
    dataPrep <- TCGAanalyze_Preprocessing(object = dataPrep, 
                                        cor.cut = 0.6,
                                        datatype = "HTSeq - Counts")                      
    dataPrep = data
    # load(paste0("./expression_HTE/DEA/TCGA_", cancer,"_HTSeq_Counts.rda"))
    dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                        geneInfo = geneInfoHT,
                                        method = "gcContent") 

    #boxplot(dataPrep, outline = FALSE)
    #boxplot(dataNorm, outline = FALSE)

    dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                    method = "quantile", 
                                    qnt.cut =  0.25)   

# https://github.com/BioinformaticsFMRP/TCGAbiolinks/issues/293
    # dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,dataSmTP], ##########TP = tumor, NT = normal
    #                             mat2 = dataFilt[,dataSmNT], ##############################  WRONG! FLIPPED
    #                             Cond1type = "Normal",
    #                             Cond2type = "Tumor",
    #                             fdr.cut = 0.01 ,
    #                             logFC.cut = 1,
    #                             method = "glmLRT")  
# Cond1type = control, Cond2type = case
    dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,dataSmNT], ##########TP = tumor, NT = normal
                                mat2 = dataFilt[,dataSmTP], ##############################  WRONG! FLIPPED
                                Cond1type = "Normal",
                                Cond2type = "Tumor",
                                fdr.cut = 0.01 ,
                                logFC.cut = 1,
                                method = "glmLRT")  

    # DEG = rownames(dataDEGs)
    # tmp = getBM(attributes = "hgnc_symbol", filters = "ensembl_gene_id", values = DEG, mart = ensembl) # 4957 out of 4979 were mapped to gene names
    # write.csv(tmp$hgnc_symbol, paste0("./DEA/", cancer, "_DEG.csv"))
    write.csv(dataDEGs, paste0("./expression_HTE/tables/", cancer, "_DEGtable.csv"))
    print(paste0(cancer, " result saved."))

}


# library(tidyr)

# # Visualization

# both = c("mat_norm", "mat_tumor")

# for (mat in both) {

#     mat <- mat[,grep("-01.-", colnames(mat))]

#     #Scale FPKM value by patients
#     mat <- mat/apply(mat,2,max)
#     mat <- t(mat)
#     mat <- cbind(separate(as.data.frame(rownames(mat)), "rownames(mat)",
#                             c(NA, "TSS", "patient", NA, "portion", "plate", "center"),  # skip var with NAs
#                             sep = "-"), mat)
#     mat$bcr <- rownames(mat)

#     # Select only one aliquot
#     mat <- as.data.frame(mat %>% dplyr::group_by(patient) %>% dplyr::slice(1))
#     tmp = strsplit(mat$bcr, "-")
#     tmp = unlist(lapply(tmp, function(x) paste(x[[1]], x[[2]], x[[3]], sep = "-")))
#     rownames(mat) <- unlist(tmp)

# }

# normal = data.frame(group = "normal", value = ptp_norm)
# tumor = data.frame(group = "tumor", value = ptp_tumor)
# ptp = rbind(normal, tumor)
# ptp$value = as.numeric(as.character(ptp$value))
# pdf("hnsc_.pdf")
# ggplot(ptp, aes(x=group, y=value, fill=group)) +  # This is the plot function
#   geom_boxplot()
#   dev.off()

