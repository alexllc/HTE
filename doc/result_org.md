# Results for publication
2020-01-18

## Main manuscript results
1. expression HTE results  `alex_lau@kwthcs-T630:~/proj/HTE/exp/2019-12-19_TCGA-pancacer_DEA_opposite_tx_dirct_exp_HTE/res/THCA`
   - Binary tx, tx direction opposite to DEA
   - FPKM-UQ normalized, linear scaled per patient
   - KM model OS imputation as outcome
 2. `#FIXME` BRCA external validation with METABRIC `exp/2020-08-21_TCGA-BRCA_METABRIC_SHC_HTE/res`
   - SHC-like validation of forests
   - binary: gene expression as treatment
   - tx direction indicated by DEA -> TCGA, CNA -> METABRIC
   - NNMIS imputed OS outcome for METABRIC and TCGA
   - no varimp overlap


## Supplementary results
1. METABRIC expression HTE `exp/2020-06-24_METABRIC_microarray_exp_HTE_NNMIS_outcome/res_with_permutate`
   - Binary tx, tx direction according to CNA
      - different nature VS DEA but highly correlated
      - there are overlaps b/t CNA and DEA but most of them are not > 50%
      - should have used both DEA and CNA to assure physiological significance
   - log$_2$ transformed array intensities as covar. Details: [METABRIC-Supplementary](https://static-content.springer.com/esm/art%3A10.1038%2Fnature10983/MediaObjects/41586_2012_BFnature10983_MOESM265_ESM.pdf)
   - NNMIS imputed outcome
2. Demographic comparison between TCGA and METABRIC cohorts `exp/2020-12-11_tx_gruop_size_test/doc/2021-01-02_demographic.md`
3. Mixed mutation and expression HTE `exp/2020-10-07_TCGA_mixed_mut_tx_exp_covar_HTE/res/2020-11-30_perm_all`
   - binary gene mutation 0/1 as tx status
   - only genes with mutation > 1% freq is considered
   - FPKM-UQ unscaled as covar
   - KM model OS imputation as outcome
4. PROPS pathway HTE `exp/2020-08-31_TCGA_pancancer_pathway_PROPS_HTE`
   - Each KEGG pathways were used in turns as treatment/covariates (i.e. when pathway A is used as tx, it's removed from the covar vector etc.)
   - NNMIS imputed outcome
   - `#FIXME` binary: only LQ as treatment group
5. Pancancer censoring proportion `exp/2020-10-07_TCGA_mixed_mut_tx_exp_covar_HTE/res`
6. BRCA external validation with METABRIC mutation only `exp/2020-06-29_TCGA-BRCA_METABRIC_SHC_only_with_varimp_validation_mut_HTE`
   - mutation as both treatment and covar
   - with varimp overlap
7. `# FIXME`Pancancer DEA indicated direction
   - only two cancer types done: BLCA, BRCA
8. Treatment group size test `exp/2020-12-11_tx_gruop_size_test/tmp/interactive_out`
   - Justification of selecting 25% as treatment or motivation to change this setting
9. `# TODO` ALDEx2 results of DE genes