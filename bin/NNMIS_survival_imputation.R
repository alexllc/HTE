# Function for imputin survival time with NNMIS on TCGA patients
# For TCGA clinical dataset, input must be a dataframe with: "bcr_patient_barcode", "age_at_initial_pathologic_diagnosis", "gender", "OS.time", "OS", "ajcc_pathologic_tumor_stage", "tumor_status" columns

library(NNMIS)
library(dplyr)

impute_with_NNMIS <- function(clin_df, type = "TCGA", only_export_obj = FALSE) {
        labels = c("[Discrepancy]","[Not Applicable]","[Not Available]","[Unknown]")
        clin_df$ajcc_pathologic_tumor_stage[which(clin_df$ajcc_pathologic_tumor_stage %in% labels)] = NA
        clin_df$type = NULL
        # Conver all variables into numeric
        for (c in colnames(clin_df)) {
            if (!is.numeric(clin_df[,c]) && c != "donorId") {
                # Convert empty strings into NAs first
                which.one <- which( levels(clin_df[,c]) == "")
                levels(clin_df[,c])[which.one] <- NA
                clin_df[,c] = sapply(sapply(clin_df[,c], as.factor), as.numeric) 
                print(paste0(c, " is converted to numeric.")) 
            }
        }
        if(all(is.na(clin_df$ajcc_pathologic_tumor_stage))) {
            message("All AJCC stage entries are NA.")
            clin_df$ajcc_pathologic_tumor_stage = NULL
            clin_df$tumor_status[floor(dim(clin_df)[1]/2)] = NA # manually removing one data point or else NNMIS will not permit using this as the auxillary variable
            tcga_imp = NNMIS(clin_df$tumor_status, 
                            xa = clin_df$age_at_initial_pathologic_diagnosis, 
                            xb = clin_df$age_at_initial_pathologic_diagnosis, 
                            time = clin_df$OS.time, 
                            event = clin_df$OS, 
                            imputeCT = T, 
                            Seed = 2020, 
                            mc.cores = 60)
            tcga_imp_surv = tcga_imp$dat.T.NNMI %>% mutate(mean = rowMeans(.))
            tcga_imp_covar = as.data.frame(sapply(tcga_imp$dat.NNMI, as.numeric)) %>% mutate(mean = rowMeans(.))
            clin_df$outcome = tcga_imp_surv$mean
            clin_df$tumor_status = tcga_imp_covar$mean
        } else {
            # Removing tumor status as an additional covariate to keep as many complete cases as possible

            # Incase there are tumor types with full records we manually remove one data point or NNMIS WILL NOT RUN
            clin_df$ajcc_pathologic_tumor_stage[floor(dim(clin_df)[1]/2)] = NA

            tcga_imp = NNMIS(clin_df$ajcc_pathologic_tumor_stage, 
                            xa = clin_df$age_at_initial_pathologic_diagnosis, 
                            xb = clin_df$age_at_initial_pathologic_diagnosis, 
                            time = clin_df$OS.time, 
                            event = clin_df$OS, 
                            imputeCT = T, 
                            Seed = 2020, 
                            mc.cores = 60)
            tcga_imp_surv = tcga_imp$dat.T.NNMI %>% mutate(mean = rowMeans(.))
            tcga_imp_covar = as.data.frame(sapply(tcga_imp$dat.NNMI, as.numeric)) %>% mutate(mean = rowMeans(.))
            clin_df$outcome = tcga_imp_surv$mean
            clin_df$ajcc_pathologic_tumor_stage = tcga_imp_covar$mean
        }
    
    if (only_export_obj) {
        return(tcga_imp)
    } else {
        return(clin_df)
    }
}