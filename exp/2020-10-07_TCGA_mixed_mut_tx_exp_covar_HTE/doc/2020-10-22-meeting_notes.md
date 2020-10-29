# Meeting notes with Larry

## Obj: 
1. Find sig genes (for prof email)
2. Immune proportion
3. where do we stop



## Dis:
1. [Larry] Try with the new mixed mut/exp HTE results (pause) we need a more systematic way of selecting genes into xg-boost
2. [Larry] Try new reference matrix (based on CPM) -> EPIC (i.e. re-run all immune proportion analysis)
    - but because EPIC is designed for TPM and the new pure immune proportion matrix is in CPM we need to use a different package
    - we can use RSEM from GDC firehose for the convoluted matrix
    - potential packages: DWLS / xCIBERSORT
    - EPIC has proportion for tumor/non-immune cell types, but for  DWLS you'd need to input the pure tumor population
    - CIBERSORTx supports absolute mode, meaning that we don't have to know the pure tumor population. 
    - Using LM22, and we don't have to generate a new convoluted matrix.

3. [Larry] LASSO / NLP (do we do this validation?)

## To-do: (for HNSC as a starter)
- [x] LM22 (from larry), LM22 is TPM based. 
    Already uploaded:
    * Server 234: `/home/yujia/Project/Immune_Deconvolution`
    * `LM22.txt`: expression matrix with pre-selected gene.
    * `LM22_source_GEPs.txt`: expression matrix with whole gene set.
- [x] scRNASeq profile for tumor types (HNSC), this can be achieved by CIBERSORTx. CIBERSORTx provides a function to generate signature file based on single cell profile.
    Already uploaded:
    * Server 234: `/home/yujia/Project/Immune_Deconvolution/HNSC_ref.txt`
    * Issue: Less cell types compared with LM22.
- [x] Run CIBERSORTx using docker
  - [x] Apply token. 
  - [x] Use results from paper "The Immune Landscape of Cancer"