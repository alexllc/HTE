## Today's objective:

1. Set up data fetch functions to get various clinical endpoints
2. Finish script for using mut as tx and exp as covar
3. compare results using different endpoints


## Procedures:

### 1. Setting up fetch functions
   1. Updated files in `/bin/TCGA_data_functions.R`, where outcome can be set as a parameter now
      1. `outParam` can be set in `fetch_clinical_data` and `impute_with_NNMIS`
      2. Can't specify download directory for `GDCdownload(m_query)`, so the working directory will be changed for these few lines only
      3. `addBatch` (for generating extra batch related variables) and `scale` (for scaling expression levels against file maximum for each expression file) can be toggled on or off in case other batch effect adjustment methods are required. `scale  = FALSE` is recommended if you wish to plot a tree, as scaling makes all splitting thresholds ~ 0.
   2. Saving files as compressed data frames are better. Whenever you want to save data frames or large tables, you should always compress it, since our servers are always full. See this [Stackexchange post](https://stackoverflow.com/questions/44477709/r-compressively-saving-large-data-frames-to-hard-drive) to see which compression format is the best. Honestly, 
   
        33404828 Oct 19 12:06 MAF_wide_BRCA.csv
        321184 Oct 19 12:11 MAF_wide_BRCA.csv.xz 
        321184 Oct 19 11:59 MAF_wide_BRCA.rds.xz 
    
    csv and rds doesn't make much of a difference, so using csv can enhance file sharing.

    [CRAN benchmarking](https://cran.r-project.org/web/packages/brotli/vignettes/benchmarks.html) seems to say that gzip is the fastest.

        279468 Oct 19 14:36 MAF_wide_BRCA.csv.gz

    Indeed, gzip compresses the most and is the fastest to load.


### 2. mut_tx_exp_covar_HTE.R script
`exp/2020-10-07_TCGA_mixed_mut_tx_exp_covar_HTE/run/mut_tx_exp_covar_HTE.R` still in progress

   1. all `./bin/*` files are automatically loaded in a loop
   2. 4 types of endpoints are processed for comparison


### 3. Endpoint comparison 

$\tau$ predictions directly?
Full HTE analysis? (with SHC and permutation?)