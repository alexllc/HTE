
# Experiment to determine whether treatment group size affects prediction results

## Motivation

We set >1% to be the mutation treatment gene requirement, but we did not get any significant (Simes P, test_calibration::diff.pval, permutation) results. So I wanted to test if we would observe such issues for genes that have already proven to be significant with our expression HTE (via permutation) would be affected by different proportions of treatment group assignment. The hypothesis is that by setting the expression level threshold so high that by setting <1% of the people get assigned into treatment group, we would loose the heterogeneity that we once observed with the 25% treatment group.

## Procedure

1. Set up dataframes of clinical and expression
2. Import list of genes to be tested (expression significant or mutated genes) <- we'll do expression first
3. Determine list of threshold to be tested
4. Loop of each threshold `thres`
    - Loop of each treatment gene
        - split patient group into treatment and control based on `thres`
        - perform permutation and test_calibration for all
        - output: Simes p, permutation, test_calibration

| gene | Simes p | permutation p vals | test_calib pvals |
| ---- | ------- | ------------------ | ---------------- |
| CDH1 | 0.006   | 0.38               | 0.0001, 0.99     |
| ...  | ...     | ...                | ...              |
- Simes p value for each of the columns above
- output table

| proportion | Simes p | permutation p vals | test_calib pvals |
| ---------- | ------- | ------------------ | ---------------- |
| 1%         | ...     | ...                | ...              |
| 2%         |         |                    |                  |
| 5%         |         |                    |                  |
| ...        |         |                    |                  |

