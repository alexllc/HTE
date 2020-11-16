
# Performing external validation with METABRIC data

## Today's objective

1. Evaluate the best endpoints for BRCA
2. Re-perform external validation with METABRIC


## Procedures

### 1. Best endpoints for BRCA

- We follow suggestions based on TCGA-CDR, so we're down to either PFI or DFI
- To compare with METABRIC data, we have to use OS regardless.
- There are many methods to make 'better' survival prediction, but for now we will just stick to the original method.
- Outcome using length **in MONTHS** without log scaling (for easier interpretations and better separation in HTE analysis)

### 2. METABRIC validation

This is going to take sometime, we have done this before but not with the mixed data. Main data preparation file is at `~/bin/bin/METABRIC_external_validation_functions.R`

Data sources: cBioPortal for the expression median z score file for both METABRIC and TCGA. location: `~/dat/METAB-TCGA-expression/exp_median_z/`

We use the already available z-score df by cBioPortal here. According to the [cBio documentation](https://docs.cbioportal.org/5.1-data-loading/data-loading/file-formats/z-score-normalization-script), z-score are calculated as:

 > Distribution based on diploid samples only: The expression distribution for unaltered copies of the gene is estimated by calculating the mean and variance of the expression values for samples in which the gene is diploid (i.e. value is "0" as reported by discrete CNA data). We call this the unaltered distribution. If the gene has no diploid samples, then its normalized expression is reported as NA.

 > Distribution based on all samples: The expression distribution of the gene is estimated by calculating the mean and variance of all samples with expression values (excludes zero's and non-numeric values like NA, Null or NaN). If the gene has samples whose expression values are all zeros or non-numeric, then its normalized expression is reported as NA.

Otherwise for every sample, the gene's normalized expression for both the profiles is reported as

$$\frac{(r-\mu)}{\sigma}$$

where $r$ is the raw expression value, and $\mu$ and $\sigma$ are the mean and standard deviation of the base population, respectively.

Since the downloaded METABRIC file did not state explicitly, also since
> cBioPortal expects z-score normalization to take place per gene

we assume the values present in the z-score file `data_mRNA_median_Zscores.txt` of METABRIC and `data_RNA_Seq_v2_mRNA_median_Zscores.txt` of TCGA are derived with this manner. Even though the formula is fairly simple, we will still use the downloaded error to prevent more error sources.

#### Reporting and best practice

According to [collins2014](https://doi.org/10.1186/1471-2288-14-40), you should include as much information about your external validation dataset as possible, so we should make some sort of table for comparison.

| features                 | TCGA | METABRIC |
| ------------------------ | ---- | -------- |
| OS/tx IQR                |      |          |
| OS/tx SD                 |      |          |
| participant demographics |      |          |
| sample  size             |      |          |
| range of predictors      |      |          |
| confidence interval      |      |          |

#### Performance measures

1. Calibration (Calibration plot, Hosmer-Lemeshow test)
2. Discrimination (c-index)
3. ROC curve (with labeled points for specificity and sensitive calculation)
4. Overall performance (Bier, $R^2$)
5. Clinical utility (decision curve analysis)

#### Action plan for external validation

- Impute missing clinical data for all studies

    - standardize imputation method (talk with Kai) -> usd the GRF built in

- Statistics to demonstrate compatibilities between internal and external datasets

    - KS test per gene to see if the underlying distribution between the two cohorts of the **same** gene are different
    - Ranges and categories of clinical, actually any continuous predictors (make sure demographic compatibility is achieved)

- Performance measures

    - Vales of R-loss criterion ("debiased.error" column)
    - Calibration similar to the built in `test_calibration` function
        - whether the internal tau prediction is well calibrated with the external tau prediction
    
#### Considerations
- Perhaps we should not perform prediction on datasets other than