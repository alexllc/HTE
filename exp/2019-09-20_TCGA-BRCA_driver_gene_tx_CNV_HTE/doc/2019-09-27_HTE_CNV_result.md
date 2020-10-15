27/09/2019 12:55

sep 23: ran 500 driver genes as treatment in turns with the other 499 genes as covariates

sep26: ran the 53 driver genes that made it to permutation (i.e. simes <0.2) as treatment, with the rest of the genes (~20,000) as covariates

To check for the duplicating issue, look out for:
1. Are all the split_half_testing result showing the same pvalues?
2. Are the variance and tau.risk of the permuted result the same? 
