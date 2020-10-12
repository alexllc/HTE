# Hybrid HTE analysis using mutation as treatment and expression and covaraites

(I guess this is the only way you can incoorperate mutationt data to it.)



For each gene, mutated = treated, un-mutated patients = contorl.
covariates are differentially expressed genes (covariate matrix can't be too big or it wouldn't know how to select appropriate features)

Should we do log transformation for random forests?
https://stats.stackexchange.com/questions/447863/log-transforming-target-var-for-training-a-random-forest-regressor

Should we report other metrics for HTE?
https://stats.stackexchange.com/questions/283760/is-cross-validation-unnecessary-for-random-forest


Examples of possible metrics:
https://www.analyticsvidhya.com/blog/2019/08/11-important-model-evaluation-error-metrics/


Slide 13 of:
https://math.usu.edu/jrstevens/bioinf/8.Forests.pdf
Presentation of heterogeneity with heatmaps

***
2020-10-08

Proceeding without selecting for variables?

***
2020-10-09

Variable selection was completed yesterday, out of 50k genes we have selected 14 using varSelRF, it's a pretty specific selection and I want to look into the details of why it's so specific.

Not all genes that I have tried using this covariate selection is heterogeneous (TP53 isn't). But some genes can be very hetereogeneous (PTEN).


***
09/10/2020 10:05

## Today's objective:

1. How to evaluate collinearity for expression data
2. How to perform parallel linear model computation (for 1.)
3. Look into details of the varSelRF algorithms
4. (if time allows) varSelRF of mutation covariates

### 1. Multiple collinearity
[Dormann review paper](https://onlinelibrary.wiley.com/doi/10.1111/j.1600-0587.2012.07348.x)

See zotero notes

### 2. Collinearity of high dimensional data (expressions)

Export mutation covariate matrix for varsel:

    "covariates_n_gp_for_varselRF.Rdata"

Perhpas that will give more 'informative' genes.

Exported dataset for most representative tree plot (Larry's script):

    "PTEN_varselRF_mut_tx_exp_covar_tumr_stat_outcome.RData"

Exported dataset for mutation

Original varselRF script cloned to current directory.


Implementation in varSelRF is too slow, consider using VSURF with using the "ranger" option to use parallel computing.

See [Genuer's guide](https://www.user2019.fr/static/pres/t25593) for full comparison of various R implementation of forests. For our use, using ranger is better than Rborist, since we have a very high variable dimension (n) but not a lot of observation (m).

    > varsel_ranger = VSURF(X, Y, RFimplem = "ranger", parallel = TRUE, verbose = TRUE)
    Thresholding step
    Growing trees.. Progress: 0%. Estimated remaining time: 6 hours, 6 minutes, 7 seconds.

We will need to do this by Rscript... Check
    
    21713

    Growing trees.. Progress: 0%. Estimated remaining time: 8 hours, 52 minutes, 48 seconds.
    Growing trees.. Progress: 0%. Estimated remaining time: 5 hours, 13 minutes, 11 seconds.
    Growing trees.. Progress: 1%. Estimated remaining time: 4 hours, 12 minutes, 11 seconds.
    Growing trees.. Progress: 1%. Estimated remaining time: 3 hours, 52 minutes, 51 seconds.

