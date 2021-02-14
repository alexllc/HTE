#!/home/alex/R-4.0.2/bin Rscript
library(optparse)
library(org.Hs.eg.db) # translate
library(dplyr)

## Run from base directory
setwd("~/project/HTE/")

## source bin scripts
bin_ls = list.files("./bin")
for (bin in bin_ls){
    source(paste0("./bin/", bin))
}

option_list = list(
    make_option(c("-p", "--pathls"), type="character", default=NULL, 
              help="comma-separated directory names containing result files", metavar="character"),
    make_option(c("-t", "--datatype", type = "character", default="expression", 
              help="covaraite data type"))
);  
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$path)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (file paths).n", call.=FALSE)
}

message("the following is your paths")
path_list <- unlist(strsplit(opt$pathls, ","))
print(path_list)

get_gene <- function(file_name, pos = NULL) {
    gene_names <- strsplit(file_name, "_|\\.")
    gene_names <- unlist(lapply(gene_names, function(x) {x[pos]}))
    return(gene_names)
}

get_partial_simes <- function(pvec, n, u_ls = NULL) {
  if (is.null (u_ls)) stop("u level value/vector must be provided.")
  psimes <- NULL
  for (u in u_ls) {
    psimes <- c(psimes, simes.partial(floor(n * u), pvec))
  } 
  return(psimes)
}

psimes_u_ls <- c(0.05, 0.1, 0.2, 0.5, 0.8)

for (exp in path_list) {
  cancerls <- list.dirs(path = paste0(exp, "/res/"), full.names = FALSE, recursive = FALSE)
  for (cancer in cancerls) {
    # permutation
    perm <- read.csv(paste0(exp, "/res/", cancer, "/", cancer, "_", opt$datatype, "_permutate_testing_result.csv"))
    # correlation
    corr <- read.csv(paste0(exp, "/res/", cancer, "/", cancer, "_", opt$datatype, "_correlation_test_result.csv"))

    #************************tau values across genes************************

    tau_files <- list.files(path = paste0(exp, "/res/", cancer, "/"), pattern = "_tau_")
    tau_tbl <- data.frame()
    for (gene_tau in tau_files) {
      tau_values <- read.csv(paste0(exp, "/res/", cancer, "/", gene_tau))
      # summary stats
      tau_sum <- summary(tau_values$tau.val)
      tau_sum_names <- names(tau_sum)
      # interpretation, CDR unit in days, present ratio/%
      max_per <- (exp(max(tau_values$tau.val)) - 1) * 100
      min_per <- (exp(min(tau_values$tau.val)) - 1) * 100
      overall <- sum(tau_values$tau.p.adjust < 0.05) / dim(tau_values)[1]

      # Calculate partial simes for patient-wise p values
      tau_psimes <- get_partial_simes(pvec = tau_values$tau.pval, n = dim(tau_values)[1], u_ls = psimes_u_ls)
      tau_sum <- c(tau_sum, max_per, min_per, overall, tau_psimes)
      names(tau_sum) <- c(tau_sum_names, "max%", "min%", "%sig_adj_tau_p", paste0("partial_simes_pval_", gsub("0\\.", "u.", psimes_u_ls)))
      tau_tbl <- rbind(tau_tbl, tau_sum)
    }
    names(tau_tbl) <- names(tau_sum)
    ## extract gene names
    gene_names <- get_gene(tau_files, pos = 3)
    tau_tbl <- cbind(gene_names, tau_tbl)

    ## translate ENSEMBL ID to hugo symbols
    translate <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys=gene_names,columns=c("SYMBOL"), keytype="ENSEMBL"))

    symbol_unique <- NULL
    for (gene in unique(translate$ENSEMBL)) {
      symbol_ls <- translate[translate$ENSEMBL == gene,]
      symbol_ls <- paste(symbol_ls$SYMBOL, collapse = "|")
      symbol_unique <- c(symbol_unique, symbol_ls)
    }

    tau_tbl <- cbind(symbol_unique, tau_tbl)

    #********************split half significance recalculation********************
    # 1. append to tau tables (sig only)
    # 2. export to separate sheet (non-sig only)

    sh_files <- list.files(path = paste0(exp, "/res/", cancer, "/"), pattern = "_split_half")
    sh_files <- sapply(tau_tbl$gene_names, function(x) sh_files[grep(x, sh_files)])
    gene_names <- get_gene(sh_files, pos = 2)
    # TODO
    corr_psimes <- data.frame()
    corr_psimes_names <- NULL
    for (file in sh_files) {
      sh_values <- read.csv(paste0(exp, "/res/", cancer, "/", file))

      # recalculate simes partial for replicates' tau values
      # DO NOT repor tau simes or partial simes p values for individual splits
      sh_psimes_ls <- NULL
      for (pcol in c("pearson.pvalue", "kendall.pvalue", "spearman.pvalue", "fisher.pval")) {
        sh_psimes <- get_partial_simes(sh_values[[pcol]], n = dim(sh_values)[1], u = psimes_u_ls)
        sh_psimes_ls <- c(sh_psimes_ls, sh_psimes)
      }
      corr_psimes <- rbind(corr_psimes, sh_psimes_ls)
    }
    corr_psimes <- cbind(gene_names, corr_psimes)
    corr_psimes_names <- c("gene", unlist(lapply(c("pearson.pvalue", "kendall.pvalue", "spearman.pvalue", "fisher.pval"), function(x) paste0(x, "_partial_simes_", gsub("0\\.", "u.", psimes_u_ls)))))
    colnames(corr_psimes) <- corr_psimes_names

    #******************************** varimp ****************************************
    # varimp enrichments
    varimp_files <- list.files(path = paste0(exp, "/res/", cancer, "/"), pattern = paste0(cancer, "_varimp_"))

    # init
    top_varimp <- data.frame()
    for (file in varimp_files) {
      varimp <- read.csv(paste0(exp, "/res/", cancer, "/", file))
      important <- varimp$varImp[varimp$varImp > quantile(varimp$varImp, 0.75)]
      names(important) <- varimp$variable[varimp$varImp > quantile(varimp$varImp, 0.75)]
      important_gene <- important[grep("ENSG*", names(important))] # only keep genes
      tx_gene_name <- get_gene(file, pos = 3)


      # check for KEGG pathway enrichments
      # Convert ensemblIDs to ENTERZ gene ID for the query
      suppressMessages(ensembl2ID <- AnnotationDbi::select(org.Hs.eg.db, keys=names(important),columns=c("ENTREZID"), keytype="ENSEMBL"))
      # pathway enriched in high varimp covaraites
      kk <- enrichKEGG(gene = ensembl2ID$ENTREZID,
                       organism = 'hsa',
                       pvalueCutoff = 0.05)
      kk <- filter(as.data.frame(kk), p.adjust < 0.05)
      k_id <- paste(kk$ID, collapse = "|")
      k_names <- paste(kk$Description, collapse = "|")
      top_varimp <- rbind(top_varimp, c(tx_gene_name, k_id, k_names))
    }
    colnames(top_varimp) <- c("gene", "KEGG_ID", "KEGG_description")

    #**************************** bind and export *****************************
    # bind: tau table, perm risks, SH corr recalc, varimp enrichment
    tau_tbl <- inner_join(tau_tbl, perm, by = c("gene_names" = "gene"))
    tau_tbl <- inner_join(tau_tbl, corr_psimes, by = c("gene_names" = "gene"))
    tau_tbl <- inner_join(tau_tbl, top_varimp, by = c("gene_names" = "gene"))

    save_name <- strsplit(exp, "/")[[1]][7]
    write.csv(tau_tbl, file = paste0("/home/alex/project/HTE/exp/2021-02-03_proj_sum/res/", save_name, "_", cancer, "_summary.csv"), row.names = FALSE)
  }
}

