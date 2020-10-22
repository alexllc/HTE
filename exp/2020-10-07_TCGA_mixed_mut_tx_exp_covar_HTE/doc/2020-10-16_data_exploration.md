# Objective: Data exploration for running hybrid HTE

## Outcome measures

Available outcome types according to TCGA-CDR publication:

Using OS or disease-specific survival (DSS) demands longer follow-up times; thus, in many clinical trials, disease-free interval (DFI) or progression-free interval (PFI) are used (Hudis et al., 2007, Punt et al., 2007; https://wiki.nci.nih.gov/plugins/servlet/mobile#content/view/24279961)

DSS = 
In cancer, the length of time after primary treatment for a cancer ends that the patient survives without any signs or symptoms of that cancer. [NIH dictionary](https://www.cancer.gov/publications/dictionaries/cancer-terms/def/disease-free-survival)

DFI = 
In cancer, the length of time after primary treatment for a cancer ends that the patient survives without any signs or symptoms of that cancer. In a clinical trial, measuring the disease-free survival is one way to see how well a new treatment works. Also called DFS, relapse-free survival, and RFS. [NIH dictionary](https://www.cancer.gov/publications/dictionaries/cancer-terms/def/disease-free-survival)


PFI =  
The length of time during and after the treatment of a disease, such as cancer, that a patient lives with the disease but it does not get worse. In a clinical trial, measuring the PFS is one way to see how well a new treatment works. Also called progression-free survival. [NIH dictionary](https://www.cancer.gov/publications/dictionaries/cancer-terms/def/pfs)

OS = 
The length of time from either the date of diagnosis or the start of treatment for a disease, such as cancer, that patients diagnosed with the disease are still alive. In a clinical trial, measuring the OS is one way to see how well a new treatment works. Also called overall survival. [NIH dictionary](https://www.cancer.gov/publications/dictionaries/cancer-terms/expand/O)

## Procedures:
1. See `exp/2020-10-07_TCGA_mixed_mut_tx_exp_covar_HTE/run/pancancer_censoring.R` for proportion exploration
2. See notes under TCGA-CDR `liu2018` in Zotero lib