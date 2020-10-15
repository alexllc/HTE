library(readxl)
library(dplyr)
library(tidyr)

# Load clinical information
cdr = read_excel("~/project/HTE/raw/TCGA-CDR-SupplementalTableS1.xlsx")
cancer_ls = unique(cdr$type)

# Loop through amount of censoring

output = data.frame(cancer_type = character(), proportion = double())

for(cancer in cancer_ls){
    cdr_sub = filter(cdr, type == cancer)
    tbl = cdr_sub %>% group_by(OS) %>% tally()
    proportion = tbl[1,2]/(tbl[1,2]+tbl[2,2])
    out_vec = cbind(cancer, proportion)
    output = rbind(output, out_vec)
}

write.csv(output, file = "../dat/proportion_of_censoring.csv", row.names = F)