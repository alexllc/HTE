#!/home/alex/R-4.0.2/bin Rscript
# The original intention of this script is to allow us to generate a report table of Tau values, correlation and permutation results easily in a standardized and efficient manner. This has been tested with: 
# expression 
# files. Since there are subtle difference in result files from all the HTE files, precaution should be taken when generating reports of other data types.

## To run this example case interactively in the console, uncomment and execute the following lines
# opt <- NULL
# opt$datatype <- "pathway" # or "expression"
# cancer <- "BRCA"
# expri <- "/home/alex/project/HTE/exp/2020-08-31_TCGA_pancancer_pathway_PROPS_HTE"

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
  library(clusterProfiler)
  KEGGdb <- download_KEGG("hsa", keggType = "pathway", keyType = "kegg")
  KEGGdb <- as.data.frame(KEGGdb$KEGGPATHID2NAME)
  KEGGdb$from  <- gsub("M", "X", KEGGdb$from)
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
  id2KEGG <- left_join(as.data.frame(pathID), KEGGdb, by = c("pathID" = "from"))
  return(unlist(translation))
}

# FIXME this will not work if too many queries are sent at once (the varimp file loop below)
# Function to retreive KEGG pathway names from KEGGREST. For some strange reasons, they won't take a list of queries even though the query size was set at 100. However, I execute the query one path at a time, I would get an error message saying "Forbidden (HTTP 403)", hence the sys.sleep break. I am guess querying more than 100 within a short period of time would not be allowed too.
kegg_all_names <- function(pathID) {
  repeats <- floor(length(pathID) / 100)
  final_rep <- length(pathID) %% 100
  pathNames <- NULL
  for (i in 1:repeats) {
    pathNames <- c(pathNames, unlist(lapply(gsub("X", "map", pathID[(1 + (i * 100) - 100):(i * 100)]), function(x) keggFind("pathway", x))))
    Sys.sleep(5)
  }
  pathNames <- c(pathNames, unlist(lapply(gsub("X", "map", pathID[(1 + (i * 100)):(i * 100 + final_rep)]), function(x) keggFind("pathway", x))))
  return(pathNames)
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
    gene_names <- get_gene(tau_files, pos = 4)
    tau_tbl <- cbind(gene_names, tau_tbl)

    ## translate ENSEMBL ID to hugo symbols
    if (opt$datatype == "pathway") {
      pathway_names <- KEGG_path_2_name(tau_tbl$gene_names)
      tau_tbl <- cbind(pathway_names, tau_tbl)
    } else {
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
    sh_files <- sapply(tau_tbl$gene_names, function(x) sh_files[grep(x, sh_files)])
    gene_names <- ifelse(opt$datatype == "pathway", get_gene(sh_files, pos = 3), get_gene(sh_files, pos = 2))

    corr_psimes <- data.frame()
    corr_psimes_names <- NULL
    for (file in sh_files) {
      sh_values <- read.csv(paste0(expri, "/res/", cancer, "/", file))

      # recalculate simes partial for replicates' tau values
      # DO NOT report tau simes or partial simes p values for individual splits
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
    # varimp enrichment (for expression HTE) and name reporting
    varimp_files <- list.files(path = paste0(expri, "/res/", cancer, "/"), pattern = paste0(cancer, "_varimp_"))
    # init
    top_varimp <- data.frame()
    for (file in varimp_files) {
      varimp <- read.csv(paste0(expri, "/res/", cancer, "/", file))
      # choosing to report only the top 25% of varimp
      important <- varimp$varImp[varimp$varImp > quantile(varimp$varImp, 0.75)]
      names(important) <- varimp$variable[varimp$varImp > quantile(varimp$varImp, 0.75)]
      message(paste0("Processing: ", which(file %in% varimp_files), " out of ", length(varimp_files)))
      if (opt$datatype == "pathway") {
        # FIXME translating mapIDs into map names with KEGGREST, will result in "Forbidden (HTTP 403)" when too many pathIDs were queries in a short period of time
        important_path <- sort(important[grep("X", names(important))], decreasing = TRUE)
        tx_path_name <- get_gene(file, pos = 4)
        top_path <- paste(unlist(lapply(strsplit(names(important_path), "\\."), function(x) x[[1]])), collapse = "|")
        # Convert pathIDs into KEGG names
        top_varimp <- rbind(top_varimp, c(tx_path_name, top_path))
        colnames(top_varimp) <- c("tx_path", "important_varimp_pathID")
        all_path <- unique(unlist(lapply(top_varimp$important_varimp_pathID, function(x) {unlist(strsplit(x, "\\|"))})))
        all_path <- cbind(all_path, kegg_all_names(all_path))
        Sys.sleep(5)
        colnames(all_path) <- c("pathID", "description")
        all_path <- as.data.frame(all_path)

        top_varimp_names <- unlist(lapply(top_varimp$important_varimp_pathID, function(x) {
          KEGG_dscp <- unlist(lapply(unlist(strsplit(x, "\\|")), function(x) {
            x = all_path$description[all_path$pathID == x]
            }))
          x <- paste(KEGG_dscp, collapse = "|")
        }))
        top_varimp <- cbind(top_varimp, top_varimp_names)
        
      } else {
        important_gene <- important[grep("ENSG*", names(important))] # only keep genes
        tx_gene_name <- ifelse(opt$datatype == "pathway", get_gene(tau_files, pos = 3), get_gene(tau_files, pos = 4))
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
        colnames(top_varimp) <- c("gene", "KEGG_ID", "KEGG_description")
      }
    }

    #**************************** bind and export *****************************
    # bind: tau table, perm risks, SH corr recalc, varimp enrichment
    tau_tbl <- inner_join(tau_tbl, perm, by = c("gene_names" = "gene"))
    tau_tbl <- inner_join(tau_tbl, corr_psimes, by = c("gene_names" = "gene"))
    tau_tbl <- inner_join(tau_tbl, top_varimp, by = c("gene_names" = "gene"))

    save_name <- strsplit(expri, "/")[[1]][7]
    write.csv(tau_tbl, file = paste0("~/project/HTE/exp/2021-02-03_proj_sum/res/", save_name, "_", cancer, "_summary.csv"), row.names = FALSE)
  }
}

