
# 03/23 Immune proportion as treatment

03/23/2020 11:51

I used EPIC to first estimate the immune cell proportion from gene expression in TPM, then the abundance of cells is used as treatment (>0.75 quantile =1, otherwise = 0). The covariates are the RPKM level of all mRNA transcripts from TCGAbiolink downloaded expression profile, along with the immune cell proportion \treatment.
