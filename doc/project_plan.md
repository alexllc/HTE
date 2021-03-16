
# Overall project TODO list

# Additional analysis
- [ ] using actual drug as treatment variables

## Parameter setting
- [ ] Test which proportion is best (with most genes passing perm test) without pre-selection
- [ ] continuous treatment
- [x] log2 transformation test: does it boost or dampen heterogeneity?
- [x] results for opposite vs in DEA treatment direction
- [x] install [warrenmcg/sleuth-ALR](https://rdrr.io/github/warrenmcg/sleuth-ALR/)

## Validation TODO
- [ ] Convert aligent using IQLR for technical concordant comparison
- [ ] re-reun expression external validation
- [ ] ext validation varimp overlap evaluation
- [ ] Gene set enrichment approaches [kloet2020](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008295)

### Compositional data
- [ ] Compare HTE results using IQLR transformed count matrix
- [ ] select covaraites/treatments based on ALDEx2 results

## Results & Presentation
- [ ] update written method section
- [ ] write scripts to generate 'reports' for results you already have

Feb 25, 2020 meeting:
- [ ] criteria for sorting visualization -> histogram etc
- [ ] genes with largest HTE -> any pathways, any famous genes, any drug targets
    - any external validationm- going on
    - add simes test (simple simes test) for tau pval
- [ ] for each tumor, prop. of genes with HTE (prop. of genes with sig.
pval_perm, pval_var, simes_test, etc)
- [ ] case study - esp. the more known genes (eg known driver genes), and/or drug targets
- [ ] pan-cancer - common and non-shared HTE genes (insight from other pan-papers?)
- [ ] Mutation results (probably we can exclude 1st)
- [ ] Ding's method mult gp -> extensionof levene test to multiple gps
- [ ] drug targets -> are any of the genes that show HTE drug targets? this may have higher clinical implications

## Additional supplementary tasks
- [ ] try X-learner with mutation data

## Advanced modifications
- PFI as outcome
- Alternative survival models (assumption scrutiny)
- Treatment continuous vs binary