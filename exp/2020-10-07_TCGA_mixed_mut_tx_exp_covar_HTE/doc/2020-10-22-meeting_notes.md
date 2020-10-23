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
    - EPIC has proportion for tumor/non-immune cell types, but for  DWLS / xCIBERSORT you'd need to input the pure tumor population

3. [Larry] LASSO / NLP (do we do this validation?)

## To-do: (for HNSC as a starter)
- [ ] LM22 (from larry)
- [ ] scRNASeq profile for tumor types (HNSC), this can be achieved by CIBERSORTx. CIBERSORTx provide a function to generate signature file based on single cell profile.
- [ ] get DWLS (R pkg)