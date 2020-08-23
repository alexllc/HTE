# This is an extension file for PROPs_score_obtain script, stored some important functions.

"obtain_tcga_data_tumor" <- function(tcga_type, download_location){
  g_query_tumor <- GDCquery(project = paste0("TCGA-", tcga_type),
                            data.category = "Transcriptome Profiling",
                            legacy = F,
                            data.type = "Gene Expression Quantification",
                            workflow.type = "HTSeq - FPKM-UQ",
                            sample.type = "Primary Tumor")
  GDCdownload(g_query_tumor, directory = paste0(download_location, "/", tcga_type, "_tumor"))
  tumor <- GDCprepare(g_query_tumor, save = FALSE, directory = paste0(download_location, "/", tcga_type, "_tumor"))
  return(tumor)
}

"obtain_tcga_data_normal" <- function(tcga_type, download_location){
  g_query_normal <- GDCquery(project =  paste0("TCGA-", tcga_type),
                            data.category = "Transcriptome Profiling",
                            legacy = F,
                            data.type = "Gene Expression Quantification",
                            workflow.type = "HTSeq - FPKM-UQ",
                            sample.type = "Solid Tissue Normal")
  if (!is.null(g_query_normal)){
    GDCdownload(g_query_normal, directory = paste0(download_location, "/", tcga_type, "_normal"))
    normal <- GDCprepare(g_query_normal, save = FALSE, directory = paste0(download_location, "/", tcga_type, "_normal"))
    return(normal)
  } else {
    return(NULL)
  }
}

"retreive_sample_type" <- function(i){ # This function is used to get the sample type code
  return(str_match(i,"^(?:[^-]*-){3}([^-(A|B|C)]*)")[2])
}

"retreieve_batch" <- function(i){ # This function is used to get the batch code
  return(str_match(i,"^(?:[^-]*-){5}(.*?)(?=-)")[2])
}

"retreieve_tss" <- function(i){ # This function is used to get the tss code
  return(str_match(i,"^(?:[^-]*-){1}([^-(A|B)]*)")[2])
}
