# This is the small script for obtaining the PROPs score for multiple TCGA-dataset

# load packages
library(optparse)
library(readr)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(sva) # (Need sva-devel version: devtools::install_github("zhangyuqing/sva-devel"))
library(PROPS)

# load needed function from function file
source("~/Network_R_project/code/PROPs_score/PROPs_score_obtain_func.r")

# Setting the args (Use optparse package)
option_list = list(
  make_option(c("-d", "--tcga_dataset_file"), type="character", default=NULL, 
              help="dataset file name (should stored as \"BRCA\\nACC\")", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="~/", 
              help="output props file directory [default= %default]", metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default="meta.txt", 
              help="output meta file name (deprecated) [default= %default]", metavar="character"),
  make_option(c("--GDCdownload"), type="character", default="~/Public_data/TCGAbiolink_download/", 
              help="GDCdownload data path [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Check whether arguments is provided
if (is.null(opt$tcga_dataset_file)){
  print_help(opt_parser)
  stop("Please provide the tcga data set id", call.=FALSE)
}

## Check whether argument is validated.
# Read the tcga_dataset_file
tcga_dataset = read_lines(opt$tcga_dataset_file)

# Check whether the dataset is valid
tcga_type_list <- c("LAML","ACC","BLCA","LGG","BRCA","CESC","CHOL","LCML","COAD","CNTL","ESCA","FPPP","GBM"
                    ,"HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSC","DLBC","MESO","MISC","OV","PAAD","PCPG","PRAD"
                    ,"READ","SARC","SKCM","STAD","TGCT","THYM","THCA","UCS","UCEC","UVM")
for (i in tcga_dataset){
  if (!(i %in% tcga_type_list)){
    print_help(opt_parser)
    stop("Please provide valid tcga data set id", call.=FALSE)
  }
}

## Main program
library(TCGAbiolinks)

# Obtain the expression matrix of selected samples
for (n in 1:length(tcga_dataset)){
  # Obtain the tcga data set
  tcgatumor <- obtain_tcga_data_tumor(tcga_dataset[n], opt$GDCdownload)
  tcganormal <- obtain_tcga_data_normal(tcga_dataset[n], opt$GDCdownload)
  
  # Summary the expression data
  exp_matrix_tumor <- SummarizedExperiment::assay(tcgatumor, "HTSeq - FPKM-UQ")
  exp_matrix_tumor <- as.data.frame(exp_matrix_tumor)
  
  # MetaLabel of tumor sample
  # tumor_label <- c(tumor_label, c(rep(paste0('Tumor-', tcga_dataset[i]), length(colnames(exp_matrix_tumor)))))
  
  if (!is.null(tcganormal)){
    # Summary the expression data
    exp_matrix_normal <- SummarizedExperiment::assay(tcganormal, "HTSeq - FPKM-UQ")
    exp_matrix_normal <- as.data.frame(exp_matrix_normal)
    
    # MetaLabel of normal sample
    # normal_label <- c(normal_label, c(rep(paste0('Normal-', tcga_dataset[i]), length(colnames(exp_matrix_normal)))))
  } else {
    exp_matrix_normal <- NULL
    print(paste0("TCGA-", tcga_dataset[n], " doesn't have normal samples, will be dropped in this analysis."))
    next
  }
  
  ######## Combine the exp matrix of different samples. (test function, will be deprecated in the future)
  # if (length(colnames(exp_matrix_normal)) != 0){
  #   exp_matrix_normal_result <- bind_cols(exp_matrix_normal_result, exp_matrix_normal)
  # }
  # 
  # exp_matrix_tumor_result <- bind_cols(exp_matrix_tumor_result, exp_matrix_tumor)
  
  ######## Applying Combat to correct the batch effect.
  # Use table() to get the count of different batch within normal and tumor samples set
  normal_batch <- as.vector(sapply(colnames(exp_matrix_normal), retreieve_batch))
  tumor_batch <- as.vector(sapply(colnames(exp_matrix_tumor), retreieve_batch))
  normal_batch_table <- table(normal_batch)
  tumor_batch_table <- table(tumor_batch)
  normal_batch_one_sample_index <- c()
  tumor_batch_one_sample_index <- c()
  for (i in 1:length(normal_batch_table)){
    if (normal_batch_table[i] == 1){
      tmp_index <- names(normal_batch_table)[i]
      normal_batch_one_sample_index <- c(normal_batch_one_sample_index, which(normal_batch %in% tmp_index))
    }
  }
  
  for (i in 1:length(tumor_batch_table)){
    if (tumor_batch_table[i] == 1){
      tmp_index <- names(tumor_batch_table)[i]
      tumor_batch_one_sample_index <- c(tumor_batch_one_sample_index, which(tumor_batch %in% tmp_index))
    }
  }
  
  # Remove the single batch.
  if (length(tumor_batch_one_sample_index) > 0){
    exp_matrix_tumor <- exp_matrix_tumor[,-tumor_batch_one_sample_index]
  }
  
  if (length(normal_batch_one_sample_index) > 0){
    exp_matrix_normal <- exp_matrix_normal[,-normal_batch_one_sample_index]
  }
  
  # Check if there is too few batches, if so, batch correction is unneeded.
  normal_singe_batch <- FALSE
  tumor_singe_batch <- FALSE
  if (length(unique(as.vector(sapply(colnames(exp_matrix_normal), retreieve_batch)))) == 1){
    normal_singe_batch <- TRUE
    print("Only one batch contained in normal set.")
  }
  
  if (length(unique(as.vector(sapply(colnames(exp_matrix_tumor), retreieve_batch)))) == 1){
    tumor_singe_batch <- TRUE
    print("Only one batch contained in tumor set.")
  }
  
  normal_count <- length(colnames(exp_matrix_normal))
  tumor_count <- length(colnames(exp_matrix_tumor))
  
  if (normal_singe_batch == TRUE | tumor_singe_batch == TRUE){ # I think I can assume that Normal samples are always in shortage.
    print("Combat will only be performed on single set.")
    # remove gene with 0 variance.
    exp_matrix_nt <- as.matrix(exp_matrix_tumor)
    var0index <- as.integer(which ((apply(exp_matrix_nt, 1, var)==0) == "TRUE"))
    exp_matrix_nt <- exp_matrix_nt[-var0index,]

    # Before perform the Combat, we need to specify a model telling combat to treat normal and tumor tissue
    # differently, otherwise we may get a inaccurtae result.
    model_tn <- data.frame("sample" = colnames(exp_matrix_nt))
    model_tn$Label <- sapply(model_tn$sample, retreive_sample_type)
    model_tn$Batch <- sapply(model_tn$sample, retreieve_batch)
    rownames(model_tn) <- model_tn$sample
    model_tn$hastumor <- model_tn$Label == "01"
    model_ComBat <- model.matrix(~hastumor, data = model_tn)
    
    exp_matrix_bc <- ComBat(dat = exp_matrix_nt, batch = model_tn$Batch, mod = model_ComBat)
    exp_matrix <- as.data.frame(exp_matrix_bc)
    
    exp_matrix$gene <- rownames(exp_matrix)
    rownames(exp_matrix) <- NULL
    exp_matrix_normal$gene <- rownames(exp_matrix_normal)
    rownames(exp_matrix_normal) <- NULL
    exp_matrix <- merge(exp_matrix_normal, exp_matrix, by = "gene", all=F, sort=F)
    rownames(exp_matrix) <- exp_matrix$gene
    exp_matrix$gene <- NULL
  } else {
    print("Combat will be performed on combined set.")
    # Combine normal tumor sample and remove gene with 0 variance.
    exp_matrix_nt <- as.matrix(cbind(exp_matrix_normal, exp_matrix_tumor))
    var0index <- as.integer(which ((apply(exp_matrix_nt, 1, var)==0) == "TRUE"))
    exp_matrix_nt <- exp_matrix_nt[-var0index,]
    
    # Before perform the Combat, we need to specify a model telling combat to treat normal and tumor tissue
    # differently, otherwise we may get a inaccurtae result.
    model_tn <- data.frame("sample" = colnames(exp_matrix_nt))
    model_tn$Label <- sapply(model_tn$sample, retreive_sample_type)
    model_tn$Batch <- sapply(model_tn$sample, retreieve_batch)
    rownames(model_tn) <- model_tn$sample
    model_tn$hastumor <- model_tn$Label == "01"
    model_ComBat <- model.matrix(~hastumor, data = model_tn)
    
    # Use ComBat function directly. (Need sva-devel version: devtools::install_github("zhangyuqing/sva-devel"))
    exp_matrix_bc <- ComBat(dat = exp_matrix_nt, batch = model_tn$Batch, mod = model_ComBat)
    exp_matrix <- as.data.frame(exp_matrix_bc)
  }
  
  ######## Use PROPS to calculate the pathway score.
  # Step 1
  # Convert the gene ID into Entrez ID. (Use bitr function to obtain a higher mapping rate)
  genes <- bitr(rownames(exp_matrix), fromType = "ENSEMBL",
                toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  # Match
  matched <-  match(rownames(exp_matrix), genes$ENSEMBL)
  exp_matrix$entrezID <-  genes$ENTREZID[matched]
  exp_matrix <- exp_matrix %>% filter(!is.na(entrezID) & entrezID != "") # Remove miss match
  dup <- duplicated(exp_matrix$entrezID) # Remove duplicate symbol row
  exp_matrix <- exp_matrix[!dup,]
  rownames(exp_matrix) <- exp_matrix$entrezID
  exp_matrix$entrezID <- NULL
  
  # Step 2
  # PROPS algorithm
  print("Start PROPS algorithm.")
  props_features_tumor <- props(as.data.frame(t(exp_matrix[,1:normal_count])), as.data.frame(t(exp_matrix[,(normal_count+1):length(exp_matrix)])))
  
  # Add a indicator variable to denote sample label.
  props_features_tumor <- props_features_tumor %>% tibble::column_to_rownames(., var = "pathway_ID")
  
  # Remove row contains inifity.
  props_features_tumor_rminf <- props_features_tumor[!is.infinite(rowSums(props_features_tumor)),]
  
  # Step 3
  # Write the file to a single file.
  write.table(props_features_tumor_rminf, file = paste0(opt$out, "props_pathway_score_", tcga_dataset[n], ".tsv"), sep = "\t")
  
  print(paste0("Finished TCGA-", tcga_dataset[n]))
}

print("Finished the script.")
# Metadata
# cat(normal_label, file=opt$metadata, append=TRUE)
# cat(" ", file=opt$metadata, append=TRUE) # To avoid the link between last normal sample and the first tumor sample
# cat(tumor_label, file=opt$metadata, append=TRUE)

# Save expression matrix (test function, will be deprecated in the future)
# if (length(colnames(exp_matrix_normal_result)) != 0){
#   result_matrix <- bind_cols(exp_matrix_normal_result, exp_matrix_tumor_result)
# } else {
#   result_matrix <- exp_matrix_tumor_result
# }
# 
# write.table(result_matrix, file = opt$out, sep = "\t")


