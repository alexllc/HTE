#!/home/alex/R-4.0.2/bin Rscript
library(optparse)
library(org.Hs.eg.db) # translate
library(dplyr)

option_list = list(
    make_option(c("-p", "--pathls"), type="character", default=NULL, 
              help="comma-separated directory names containing result files", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="out.csv", 
              help="output file name [default= %default]", metavar="character"),
    make_option(c("-c", "--cancerls"), type = "character", default=NULL,
              help="cancer type to process"),
    make_option(c("-d", "--direction"), type = "character", default="along", 
              help="bottom or top 25\\% used as treatment group."),
    make_option(c("-t", "--datatype", type = "character", default="expression", 
              help="covaraite data type"))
);  
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$path)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (file paths).n", call.=FALSE)
}

if (is.null(opt$cancer)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (cancer type).n", call.=FALSE)
}


message("the following is your paths")
path_list <- unlist(strsplit(opt$pathls, ","))
print(path_list)

message("the following is your cancer types")
print(opt$cancerls)

get_gene <- function(file_name, pos = NULL) {
    gene_names <- strsplit(file_name, "_|\\.")
    gene_names <- unlist(lapply(gene_names, function(x) {x[pos]}))
    return(gene_names)
}


for (cancer in cancerls) {
  
  for (exp in path_list) {
    # permutation
    perm <- read.csv(paste0(exp, "/res/", cancer, "/", cancer, "_", opt$datatype, "_permutate_testing_result.csv"))

    write.xlsx(perm, file = "summary_output.xlsx", sheetName = paste0("DE", opt$direction, "_permutation"), append = FALSE)
  
    # correlation
    corr <- read.csv(paste0(exp, "/res/", cancer, "/", cancer, "_", opt$datatype, "_correlation_test_result.csv"))
    
    write.xlsx(perm, file = "summary_output.xlsx", sheetName = paste0("DE", opt$direction, "_correlation"), append = TRUE)

    # tau values across genes
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

      tau_sum <- c(tau_sum, max_per, min_per, overall)
      names(tau_sum) <- c(tau_sum_names, "max%", "min%", "%sig_adj_tau_p")
      tau_tbl <- rbind(tau_tbl, tau_sum)

    }
    names(tau_tbl) <- names(tau_sum)

    ## extract gene names
    gene_names <- get_gene(tau_files)
    tau_tbl <- cbind(gene_names, tau_tbl)

    ## translate ENSEMBL ID to hugo symbols
    translate <- AnnotationDbi::select(org.Hs.eg.db, keys=gene_names,columns=c("SYMBOL"), keytype="ENSEMBL")

    symbol_unique <- NULL
    for (gene in unique(translate$ENSEMBL)) {
      symbol_ls <- translate[translate$ENSEMBL == gene,]
      symbol_ls <- paste(symbol_ls$SYMBOL, collapse = "|")
      symbol_unique <- c(symbol_unique, symbol_ls)
    }

    tau_tbl <- cbind(symbol_unique, tau_tbl)

    write.xlsx(tau_tbl, file = "summary_output.xlsx", sheetName = paste0("DE", opt$direction, "_tau_values"), append = TRUE)

    # split half significance recalculation
    sh_files <- list.files(path = paste0(exp, "/res/", cancer, "/"), pattern = "_split_half")
    gene_names <- get_gene(sh_files, pos = 2)

    # TODO
    for (i in 1:length(sh_files)) {
      sh_values <- read.csv(paste0(exp, "/res/", cancer, "/", sh_files[i]))
      gene_name <- gene_names[i]

      # simes re calculate as proportion significant
        for (j in 1:10) {
            sh_a <- read.csv(paste0(exp, "/res/", "/", cancer, "/", cancer, "_", gene_name, "_observation_", j, "_result_a.csv"))
            sh_b <- read.csv(paste0(exp, "/res/", "/", cancer, "/", cancer, "_", gene_name, "_observation_", j, "_result_b.csv"))
            rep_a[,i] <- sh_a$V1
            rep_b[,i] <- sh_b$V1
        }

        
      # Identify the case 2 genes

      
      # correlation


    }

    # varimp enrichments

    # pathway activation
  }

}