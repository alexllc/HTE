## Today's objective:

1. Evaluate the best endpoints for BRCA
2. Re-perform external validation with METABRIC


## Procedures:

### 1. Best endpoints for BRCA
- We follow suggestions based on TCGA-CDR, so we're down to either PFI or DFI
- To compare with METABRIC data, we have to use OS regardless.
- There are many methods to make 'better' survival prediction, but for now we will just stick to the original method.
- Outcome using length **in MONTHS** without log scaling (for easier interpretations and better separation in HTE analysis)

### 2. METABRIC validation
This is going to take sometime, we have done this before but not with the mixed data. Main data preparation file is at `~/bin/bin/METABRIC_external_validation_functions.R`

Data sources: cBioPortal for the expression median z score file for both METABRIC and TCGA. location: `~/dat/METAB-TCGA-expression/`