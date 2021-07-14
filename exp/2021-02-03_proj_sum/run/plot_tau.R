# plot tau values

## Run from base directory
setwd("~/project/HTE/")

## source bin scripts
bin_ls = list.files("./bin")
for (bin in bin_ls){
    suppressMessages(source(paste0("./bin/", bin)))
}

option_list = list(
    make_option(c("-t", "--tau"), type="character", default=NULL, 
              help="comma-separated directory names containing tau value files", metavar="character"),
    make_option(c("-", "--", type = "character", default="expression", 
              help="covaraite data type"))
);  
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);