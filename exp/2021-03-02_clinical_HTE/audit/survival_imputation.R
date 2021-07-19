# get OS time from TCGA-CDR
cdr <- read_excel("./dat/TCGA-CDR-SupplementalTableS1.xlsx")
outParam = "OS"
col_vec =  c("bcr_patient_barcode", "type", outParam, paste0(outParam, ".time"))
clinical_dat <- dplyr::select(cdr, col_vec)

# labels <- c("[Discrepancy]","[Not Applicable]","[Not Available]","[Unknown]")
# clinical_dat$ajcc_pathologic_tumor_stage[which(clinical_dat$ajcc_pathologic_tumor_stage %in% labels)] <- NA
clinical_dat$type = NULL
clinical_dat[clinical_dat == "#N/A"] <- NA # sometimes missing entries are not well formatted
clinical_dat <- subset(clinical_dat, !is.na(get(paste0(outParam, ".time"))) & !is.na(get(outParam)))
colnames(clinical_dat)[colnames(clinical_dat)=="bcr_patient_barcode"] = "donorId"
colnames(clinical_dat)[colnames(clinical_dat) == outParam | colnames(clinical_dat) == paste0(outParam, ".time")] = c("out_param", "out_param_time")
clinical_dat = as.data.frame(clinical_dat)

# get clinical data from BCR Biotab




csf <- causal_survival_forest(
                                X = X,
                                Y = wds$out_param_time,
                                W = wds$PEMETREXED,
                                D = wds$out_param)