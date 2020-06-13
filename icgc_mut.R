library(UCSCXenaTools)
library(dplyr)
library(tidyverse)
data(XenaData)
head(XenaData)
XenaGenerate(subset = XenaHostNames=="icgcHub") %>% 
  XenaFilter(filterDatasets = "SNV") %>% 
  XenaFilter(filterCohorts = "donor") %>% XenaFilter(filterDatasets = "nonUS") -> df_snv

XenaQuery(df_snv) %>% XenaDownload() -> xe_download
mut = XenaPrepare(xe_download)

XenaGenerate(subset = XenaHostNames=="icgcHub") %>% 
  XenaFilter(filterCohorts = "donor") %>% 
  XenaFilter(filterDatasets = "survival|phenotype") -> df_clinical

XenaQuery(df_clinical) %>% XenaDownload() -> xe_download
clinical = XenaPrepare(xe_download)
OS = clinical[[1]]
pheno = clinical[[2]]
  colnames(pheno)[colnames(pheno) == "_primary_site"] = "primary_site"

cancer_list = c("Blood", "Brain", "Breast", "Liver", "Kidney", "Head and neck", "Prostate", "Lung", "Stomach", "Nervous System", "Pancreas", "Colorectal", "Esophagus", "Skin", "Ovary", "Uterus", "Bladder", "Bone", "Gall Bladder", "Cervix", "Mesenchymal", "Nasopharynx")

clincol = c('sampleID', 'donor_age_at_diagnosis', 'donor_sex', 'project_code')

for (c in cancer_list) {

  df_clin <- pheno[pheno$primary_site == c,]

  # 1 Combine and tidy cancer specific datasets
  ## Coluumns we need: donorID, age, sex, (stage)

  df_clin = dplyr::select(df_clin, clincol) %>% left_join(OS, by = c("sampleID" = "icgc_donor_id")) %>% 
  # filter(str_detect(project_code, "-US$", negate = T)) %>% 
  filter(!is.na(OS.time))

  # Discard US samples, because they're incomplete records of TCGA data
  # mutsub = mut[mut$sample %in% df_clin$sampleID,]
  # mutsub = mutsub %>% filter(effect != "synonymous_variant") %>% filter(!is.na(gene))
  # mutsub = mutsub %>% group_by(sample) %>% group_by(gene) %>% add_tally() %>% dplyr::select(c("sample", "gene", "n")) %>% unique() %>% spread(gene, n, fill = 0)
  message(paste0("Number of available and valid non-US samples for ", c, ": ", length(which(df_clin$sampleID %in% mut$sample))))

}

df_comb = left_join(df_clin, nonUS, by = c("sampleID"= "sample"))