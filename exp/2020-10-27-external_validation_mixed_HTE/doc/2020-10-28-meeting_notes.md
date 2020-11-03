# Meeting with Prof So and Larry

## Topics to discuss

1. Mutation as treatment, exp as covar analysis, pval, tau
2. External validation with METABRIC and TCGA
3. Immune proportion HTE
4. Uniformity of source code and modeling parameters
5. Manuscript writing, discussion, interpretation of data


## Discussion

1. Logitp requires independent p values, simes robust to positive dependencies
   1. recalculate Simes pvalues for
   2. Perform SHC (do anyways) and permutation
    - Complete cases vs impute with GRF's updated missing value function
    - Ding's method (from Berkley), comparing variance between treated and untreated group (non-parametric method, Levene's test)
    - stability vs heterogeneity detection. SHC significant only means stability, third factor leads to correlation in SHC for treatment without real effect, any systemic factors leading to significant SHC other than true tau effects
2. Aspects of external validation: stability of models, test calibration, direct overlap
3. Use already available immune proportion


## To-do

- [ ] Perform SHC, permute for mixed HTE
- [ ] further methods to perform external validation, email prof So about your ideas
- [ ] discuss with everybody in the project to use the same set of codes.
- [ ] Impute rather than omit missing data janssen2010, cismondi2013


## Tentative procedures

### External validation

Basic requirement is 