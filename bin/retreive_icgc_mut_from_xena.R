library(UCSCXenaTools)
library(dplyr)
library(tidyverse)
data(XenaData)
head(XenaData)
library(NNMIS)

XenaGenerate(subset = XenaHostNames=="icgcHub") %>% 
  XenaFilter(filterDatasets = "SNV") %>% 
  XenaFilter(filterCohorts = "donor") %>% XenaFilter(filterDatasets = "nonUS") -> df_snv

XenaQuery(df_snv) %>% XenaDownload() -> xe_download
mut = XenaPrepare(xe_download)

XenaGenerate(subset = XenaHostNames=="icgcHub") %>% 
  XenaFilter(filterCohorts = "donor")%>%
  XenaFilter(filterDatasets = "survival|phenotype") -> df_clinical

XenaQuery(df_clinical) %>% XenaDownload() -> xe_download
clinical = XenaPrepare(xe_download)
OS = clinical[[1]]
pheno = clinical[[2]]
colnames(pheno)[colnames(pheno) == "_primary_site"] = "primary_site"

cancer_list = c("Blood", "Brain", "Breast", "Liver", "Kidney", "Head and neck", "Prostate", "Lung", "Stomach", "Nervous System", "Pancreas", "Colorectal", "Esophagus", "Skin", "Ovary", "Uterus", "Bladder", "Bone", "Gall Bladder", "Cervix", "Mesenchymal", "Nasopharynx")

clincol = c("sampleID", "donor_age_at_diagnosis", "donor_sex", "project_code", "donor_tumour_stage_at_diagnosis", "donor_tumour_stage_at_diagnosis_supplemental", "donor_tumour_staging_system_at_diagnosis")

for (c in cancer_list) {

  df_clin <- pheno[pheno$primary_site == c,]

  # 1 Combine and tidy cancer specific datasets
  ## Coluumns we need: donorID, age, sex, (stage)

  df_clin = dplyr::select(df_clin, clincol) %>% left_join(OS, by = c("sampleID" = "icgc_donor_id"))
  df_clin = df_clin %>% filter(str_detect(project_code, "-US$", negate = T)) %>% filter(!is.na(OS.time))

  # Format staing information for survival imputation
  tmp = df_clin$donor_tumour_stage_at_diagnosis[grep("^T[[:digit:]]", df_clin$donor_tumour_stage_at_diagnosis)]
  tmp = strsplit(tmp, "N")
  tmp = sapply(tmp, function(x) x[1][1])
  tmp = gsub("T", "", tmp)
  df_clin$donor_tumour_stage_at_diagnosis[grep("^T[[:digit:]]", df_clin$donor_tumour_stage_at_diagnosis)] = tmp

  df_clin$donor_tumour_stage_at_diagnosis[dim(df_clin[1]) / 2] = NA # manually create missing data
  icgc_imp = NNMIS(df_clin$donor_tumour_stage_at_diagnosis, 
                            xa = df_clin$donor_age_at_diagnosis, 
                            xb = df_clin$donor_age_at_diagnosis, 
                            time = df_clin$OS.time, 
                            event = df_clin$OS, 
                            imputeCT = T, 
                            Seed = 2020, 
                            mc.cores = 60)
  icgc_imp_surv = icgc_imp$dat.T.NNMI %>% mutate(mean = rowMeans(.))

  icgc_imp_covar = as.data.frame(sapply(icgc_imp$dat.NNMI, as.numeric)) %>% mutate(mean = rowMeans(.))
  df_clin$outcome = icgc_imp_surv$mean
  df_clin$donor_tumour_stage_at_diagnosis = icgc_imp_covar$mean

  # Discard US samples, because they're incomplete records of TCGA data 
  # removing filter(effect != "synonymous_variant") 
  mutsub = mut[mut$sample %in% df_clin$sampleID,]
  mutsub = mutsub %>% filter(!is.na(gene)) %>% group_by(sample) %>% group_by(gene) %>% add_tally() %>% dplyr::select(c("sample", "gene", "n")) %>% unique() %>% spread(gene, n, fill = 0)
  mutsub = mutsub %>% select(-(colnames(mutsub)[grep("^\\d{,2}-[[:alpha:]]{3}$", colnames(mutsub))])) # for some reason the pipe above wouldn't work

  message(paste0("Number of available and valid non-US samples for ", c, ": ", length(which(df_clin$sampleID %in% mut$sample))))

  # submut = mut[mut$sample %in% df_clin$sampleID,]
  # submut = dplyr::filter(submut, !is.na(gene)) %>% group_by(sample) %>% add_count(gene) %>% ungroup() %>% distinct()%>% spread(gene, n)

  clin_mut = left_join(df_clin, mutsub, by = c("sampleID" = "sample"))
  
  proj_name = strsplit(clin_mut$project_code[1], '-')
  proj_name = proj_name[[1]][1]
  write.csv(clin_mut, file = paste0("ICGC_", proj_name, "_df.csv"), row.names = F)

# select(-all_of(c( "chr", "start", "end", "reference", "alt", "effect", "Amino_Acid_Change"))) 
}

# df_comb = left_join(df_clin, nonUS, by = c("sampleID"= "sample"))


# for exploration

