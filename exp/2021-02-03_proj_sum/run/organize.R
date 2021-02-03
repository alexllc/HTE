#!/home/alex/R-4.0.2/bin Rscript
library("optparse")

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


for (cancer in cancerls) {
  
  for (exp in path_list) {
    # permutation
    perm <- read.csv(paste0(exp, "/res/", cancer, "/", cancer, "_", opt$datatype, "_permutate_testing_result.csv"))

    write.xlsx(perm, file = "summary_output.xlsx", sheetName = paste0("DE", opt$direction, "_permutation"), append = FALSE)
  
    # correlation
    corr <- read.csv(paste0(exp, "/res/", cancer, "/", cancer, "_", opt$datatype, "_correlation_test_result.csv"))
        write.xlsx(perm, file = "summary_output.xlsx", sheetName = paste0("DE", opt$direction, "_permutation"), append = FALSE)


    # tau values across genes

    # varimp enrichments

    # pathway activation
  }

}