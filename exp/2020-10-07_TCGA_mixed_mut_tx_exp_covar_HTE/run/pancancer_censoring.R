library(readxl)
library(dplyr)
library(tidyr)


################################################################################
# Load clinical information
cdr = read_excel("~/project/HTE/raw/TCGA-CDR-SupplementalTableS1.xlsx")
cancer_ls = unique(cdr$type)

endpt_ls = c("OS", "DSS", "DFI", "PFI")

# Loop through amount of censoring

for (endpt in endpt_ls) {
    output = data.frame(cancer_type = character(), proportion = double())

    for(cancer in cancer_ls){
        cdr_sub = filter(cdr, type == cancer)
        tbl = cdr_sub %>% group_by(get(endpt)) %>% tally()
        proportion = tbl[1,2]/(tbl[1,2]+tbl[2,2])
        out_vec = cbind(cancer, proportion)
        output = rbind(output, out_vec)
    }
    assign(paste0(endpt, "_outdput"), output)    
}
