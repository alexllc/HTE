#!/home/alex/R-4.0.2/bin Rscript
# The original intention of this script is to allow us to generate a report table of Tau values, correlation and permutation results easily in a standardized and efficient manner. This has been tested with: 
# expression 
# files. Since there are subtle difference in result files from all the HTE files, precaution should be taken when generating reports of other data types.

## To run this example case interactively in the console, uncomment and execute the following lines
# opt <- NULL
# opt$datatype <- "pathway" # or "expression"
# cancer <- "BRCA"
# expri <- "/home/alex/project/HTE/exp/2020-06-20_TCGA_pancancer_DEA_indicated_tx_dirct_.25_gp_as_tx_exp_HTE"

library(optparse)
library(org.Hs.eg.db) # translate gene ENSEMBL ID to HugoSymbols
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

if (opt$datatype == "pathway") {
  download.file("http://rest.kegg.jp/list/pathway", "path_name_ls.txt")
  KEGGdb <- read.table("path_name_ls.txt", sep = "\t")
} 

if (is.null(opt$path)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (file paths).n", call.=FALSE)
}

message("the following is your paths")
path_list <- unlist(strsplit(opt$pathls, ","))
print(path_list)

# Function to extract gene name from .csv file names
get_gene <- function(file_name, pos = NULL) {
    gene_names <- strsplit(file_name, "_|\\.")
    gene_names <- unlist(lapply(gene_names, function(x) {x[pos]}))
    return(gene_names)
}

# Recalculate a list of partial Simes values based on the simes.partial function in the HTE_validation_function.R
get_partial_simes <- function(pvec, n, u_ls = NULL) {
  if (is.null (u_ls)) stop("u level value/vector must be provided.")
  psimes <- NULL
  for (u in u_ls) {
    psimes <- c(psimes, simes.partial(floor(n * u), pvec))
  } 
  return(psimes)
}

# FIXME this downloaded DB is missing a lot of names
KEGG_path_2_name <- function(pathID = NULL) {
  pathID <- gsub("X", "path:map", pathID)
  id2KEGG <- left_join(as.data.frame(pathID), KEGGdb, by = c("pathID" = "V1"))
  return(id2KEGG$V2)
}

# determine the u levels to be used when re-calculating SH partial Simes values
psimes_u_ls <- c(0.05, 0.1, 0.2, 0.5, 0.8)

for (expri in path_list) { # 'exp' and 'expr' are reserved names in r-base
  cancerls <- list.dirs(path = paste0(expri, "/res/"), full.names = FALSE, recursive = FALSE)
  for (cancer in cancerls) {
    message(paste0("Processing cancer: ", cancer))
    # import permutation and correlation result files
    if (opt$datatype == "pathway") {
      perm <- read.csv(paste0(expri, "/res/", cancer, "/", "PROPS_", cancer, "_expression_permutate_testing_result.csv"))
      corr <- read.csv(paste0(expri, "/res/", cancer, "/", "PROPS_", cancer, "_expression_correlation_test_result.csv"))
    } else {
      perm <- read.csv(paste0(expri, "/res/", cancer, "/", cancer, "_", opt$datatype, "_permutate_testing_result.csv"))
      corr <- read.csv(paste0(expri, "/res/", cancer, "/", cancer, "_", opt$datatype, "_correlation_test_result.csv"))
    }

    #************************tau values across genes************************

    # import tau value files and summarize the tau distributions for each treatment variable
    tau_files <- list.files(path = paste0(expri, "/res/", cancer, "/"), pattern = "_tau_")
    tau_tbl <- data.frame()
    for (gene_tau in tau_files) {
      tau_values <- read.csv(paste0(expri, "/res/", cancer, "/", gene_tau))
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
    if (opt$datatype == "pathway") {
      gene_names <- get_gene(tau_files, pos = 4)
    } else { 
      gene_names <- get_gene(tau_files, pos = 3)
    }
    tau_tbl <- cbind(gene_names, tau_tbl) # bind treatment gene names

    ## translate ENSEMBL ID to hugo symbols if not already
    if (opt$datatype == "pathway") {
      pathway_names <- KEGG_path_2_name(tau_tbl$gene_names)
      tau_tbl <- cbind(pathway_names, tau_tbl)
    } else if (all(grepl("ENSG", tau_tbl$gene_names))) {
        translate <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys=gene_names,columns=c("SYMBOL"), keytype="ENSEMBL"))
        symbol_unique <- NULL
        for (gene in unique(translate$ENSEMBL)) {
          symbol_ls <- translate[translate$ENSEMBL == gene,]
          symbol_ls <- paste(symbol_ls$SYMBOL, collapse = "|")
          symbol_unique <- c(symbol_unique, symbol_ls)
        }
        tau_tbl <- cbind(symbol_unique, tau_tbl)
    }

    #********************split half significance recalculation********************
    # 1. append to tau tables (sig only)
    # 2. export to separate sheet (non-sig only)

    sh_files <- list.files(path = paste0(expri, "/res/", cancer, "/"), pattern = "_split_half")
    sh_files <- sapply(tau_tbl$gene_names, function(x) sh_files[grep(paste0("_", x, "_"), sh_files)])

    corr_psimes <- data.frame()
    corr_psimes_names <- NULL
    for (file in sh_files) {
      sh_values <- read.csv(paste0(expri, "/res/", cancer, "/", file))

      # recalculate simes partial for replicates' tau values
      # DO NOT report tau simes or partial simes p values for individual splits
      sh_psimes_ls <- NULL
      for (pcol in c("pearson.pvalue", "kendall.pvalue", "spearman.pvalue", "fisher.pval", "t.test.a.pval", "t.test.b.pval")) {
        sh_psimes <- get_partial_simes(sh_values[[pcol]], n = dim(sh_values)[1], u = psimes_u_ls)
        sh_psimes_ls <- c(sh_psimes_ls, sh_psimes)
      }
      corr_psimes <- rbind(corr_psimes, sh_psimes_ls)
    }
    if (opt$datatype == "pathway") {
      gene_names <- get_gene(sh_files, pos = 3)
    } else { 
      gene_names <- get_gene(sh_files, pos = 2)
    }
    corr_psimes <- cbind(gene_names, corr_psimes) # bind treatment gene names

    # rename all correlation Simes p value calculation
    corr_psimes_names <- c("gene", unlist(lapply(c("pearson.pvalue", "kendall.pvalue", "spearman.pvalue", "fisher.pval", "t.test.a.pval", "t.test.b.pval"), function(x) paste0(x, "_partial_simes_", gsub("0\\.", "u.", psimes_u_ls)))))
    colnames(corr_psimes) <- corr_psimes_names

    #******************************** varimp ****************************************
    # varimp enrichment (for expression HTE) and name reporting
    varimp_files <- list.files(path = paste0(expri, "/res/", cancer, "/"), pattern = paste0(cancer, "_varimp_"))
    # init
    top_varimp <- data.frame()
    for (file in varimp_files) {
      varimp <- read.csv(paste0(expri, "/res/", cancer, "/", file))
      # choosing to report only the top 25% of varimp and include translation
      important <- varimp$varImp[varimp$varImp > quantile(varimp$varImp, 0.75)]
      names(important) <- varimp$variable[varimp$varImp > quantile(varimp$varImp, 0.75)]
      if (opt$datatype == "pathway") {
        # selection of top 1/4
        important_path <- sort(important[grep("X", names(important))], decreasing = TRUE)
        tx_path_name <- get_gene(file, pos = 4)
        important_path <- unlist(lapply(strsplit(names(important_path), "\\."), function(x) x[[1]]))

        # concatenate the top 25% covar of *EACH* treatment gene into *ONE* chatacter
        top_path <- paste(important_path, collapse = "|")

        # Convert pathIDs into KEGG names
        top_varimp <- rbind(top_varimp, c(tx_path_name, top_path, paste(KEGG_path_2_name(important_path), collapse = "|")))
        colnames(top_varimp) <- c("gene", "important_varimp_pathID", "important_varimp_pathName")
        
      } else {
        # filter covarates that are *not* transcripts, they could be clinical variables
        tx_gene_name <- get_gene(file, pos = 3)

        # check for KEGG pathway enrichments
        if (all(grepl("ENSG", tau_tbl$gene_names))) {
          important_gene <- important[grep("ENSG*", names(important))] # only keep genes
          # Convert ensemblIDs to ENTERZ gene ID for the query
          suppressMessages(ensembl2ID <- AnnotationDbi::select(org.Hs.eg.db, keys=names(important),columns=c("ENTREZID"), keytype="ENSEMBL")) 
        } else {
          important_gene <- important[!grepl("_", names(important))]
          # convert Hugo gene symbol to ENTREZID
          suppressMessages(ensembl2ID <- AnnotationDbi::select(org.Hs.eg.db, keys=names(important),columns=c("ENTREZID"), keytype="SYMBOL"))

        }
        # pathway enriched in high varimp covaraites
        kk <- enrichKEGG(gene = ensembl2ID$ENTREZID,
                        organism = 'hsa',
                        pvalueCutoff = 0.05)
        kk <- filter(as.data.frame(kk), p.adjust < 0.05)
        k_id <- paste(kk$ID, collapse = "|")
        k_names <- paste(kk$Description, collapse = "|")
        top_varimp <- rbind(top_varimp, c(tx_gene_name, k_id, k_names))
        colnames(top_varimp) <- c("gene", "KEGG_ID", "KEGG_description")
      }
    }
    
    #**************************** test_calibration *****************************
    test_calib <- list.files(path = paste0(expri, "/res/", cancer, "/"), pattern = "_calibration_")
    test_calib <- read.csv(paste0(expri, "/res/", cancer, "/", test_calib))
    test_calib <- dplyr::select(test_calib, all_of(c("gene", last(colnames(test_calib), n = 4))))
    test_calib$gene <- unlist(lapply(strsplit(perm$gene, "\\."), function(x) x[1]))

    #**************************** bind and export *****************************
    # bind: tau table, perm risks, SH corr recalc, varimp enrichment
    if (opt$datatype == "pathway") 
      perm$gene <- unlist(lapply(strsplit(perm$gene, "\\."), function(x) x[1]))

      sum_tbl <- inner_join(tau_tbl, perm, by = c("gene_names" = "gene"))
      sum_tbl <- inner_join(sum_tbl, corr_psimes, by = c("gene_names" = "gene"))
      sum_tbl <- inner_join(sum_tbl, test_calib, by = c("gene_names" = "gene"))
      sum_tbl <- inner_join(sum_tbl, top_varimp, by = c("gene_names" = "gene"))

    save_name <- strsplit(expri, "/")[[1]][7]
    write.csv(sum_tbl, file = paste0("~/project/HTE/exp/2021-02-03_proj_sum/res/", save_name, "_", cancer, "_summary.csv"), row.names = FALSE)
  }
}

