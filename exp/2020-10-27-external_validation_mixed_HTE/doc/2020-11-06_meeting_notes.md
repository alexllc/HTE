
# Meetings notes with Prof So and Kai

## Items to discuss

- Unifying scripts for HTE analysis
- Survival imputation methods (using KM or other auxillary variable methods, if we're re-running should we use the survival forest instead?)
- external validation

Email except:
> Referring to our latest progress on TCGA below, we may need to come up with a prompt solution to some inconsistency of the models. Also, method for survival imputation may need to be discussed. On the other hand, since ukbb and TCGA are rather different datasets, some differences in parameter settings may be okay.

## meeting minutes

### 1. 5000 vs 2000 doesn't matter?

- COVID treatment proportion too little, but 1/4 tx is alright
- expression not too sparse
- each DEA used as tx (too much gene too computationally expensive)
- no colinear issues if they're both covariates, but if the treatment is closely associated with treatment, then it may control for the effects.
  - we may need to assess colinearity, filter correlation $R^2$ > 0.6 then exclude
  - but we don't necessarily need do filter? ***binary tx will dampen the colinearity issue*** (write in justification)
  - super small effect of each gene will not result in large colinearity
- use 5000 tree model for mutation as treatment HTE analysis
- 2000 tree is default setting

### 2. Survival imputation

- use KM

### 3. External validation

- there are flexibility, depends on your model, it's hard to be the same as other studies
- hard to do external validation if you do unsupervised similar to SHC?
- ext validation mainly for consistency and stability (bi clustering model)
- risk (we don't have standard to compare, we don't know about the variance) but comparing risks does not give you much information 
- compare external dataset predicted tau vs self-built forest 
Calibration test
$$
    lm(y~x+0)
    Z = abs(coef - 1)/SE(coeff)
$$
- reproduction of results i.e. permutated gene significant

## To-do

- [ ] Run mixed mutation, expression HTE for BRCA (use z score for all expression for consistency)
- [ ] External validation strategies finalize and run script!!
- [ ] Summary table of results
