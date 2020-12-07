
# Side notes while I wrote the script for the functions

## Random remarks

1. I think it's not right to say that the model successfully replicated results based on a simple regression between the original predicted tau values and stats with regression. However, if the models are significant (given by mean fores prediction p values), then we can say even though there isn't a significantly different heterogeneous treatment effect observed among the current patient cohort, the predicted tau values at least tell us about the random fluctuations people have on their average treatment effect.
2. So I will still perform the OOB dataset applied to causal forest validation with correlation, but together with the test calibration statistics.
3. No need to further subset the covariate matrix as it was not subset before in expression-treatment expression expression-covariate HTE analysis.


## Temp todo checklist for the function script

- [x] add cross-model consistency checking function
- [x] code and generate overall significance (compress into one value per gene per patient so that we can get better hold of what is and what isn't significant)
- [x] Add in parts where it will do what external validation have been done for the past year -> modified quantile t test to SIGN test
- [ ] Actually run the external validation script

Using `SIGN.test` from the `BDSA` package as we do not assume symmetry about median for non-parametric dependent sample test Wilcoxon or normality under T tests.
