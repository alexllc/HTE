require(methods)
require(doParallel)
require(missForest)
require(data.table)
require(grf)
require(dplyr)
require(stringr)

registerDoParallel(10)

source('causal_inference_models.R')
source('hte.validation.R')

setwd("/home/kai/data/UKBB")

handler <- function() {
    # Save debugging info to file last.dump.rda
    # the dump file can be loaded in R for debugging.
    dump.frames(to.file = TRUE)
    # Quit R with error status
    q(status = 1)
}
options(error = handler)

# if the generic function existis, i.e., append, we should use the same argument names as exist.
# otherwise, we need to create a new one.
# setGeneric("append") # this declear may be omited since append function already exists.
setMethod("append", signature(x = "data.frame", values = "vector"),
    function(x, values) {
        if(ncol(x) != length(values)){
            stop('dimension of data.frame does not match with the length of values.')
        }else if(typeof(x) == 'list'){
            x[nrow(x)+1, ] <- values
            x
        }else{
            x[nrow(x)+1, ] <- as.list(values)
            x
        }
    }
)

# Note: updated by Kai (Nov-15-2020) 
# To accommodate the inclusion of T1DM, "2443-1.0" is introduced to represent "T2DM"; 
# Meanwhile, "2443-0.0" represents "T1DM".
diabetes_exc_HbA1c <- c('2443-0.0','2443-1.0','30750-0.0')
diabetes_inc_HbA1c <- c('2443-1.0','2443-1.0')
drink_status <- c('20117-0.0', '20117-0.1', '20117-0.2', '20117-0.3')
smoking_status <- c('20116-0.0', '20116-0.1', '20116-0.2', '20116-0.3')
ethnic <- c('21000-0.0', '21000-0.1', '21000-0.2', '21000-0.3', '21000-0.4', '21000-0.5','21000-0.6')
bmi <- c('21001-0.0', '23099-0.0', '48-0.0', 'whr') # added whr (28/10/2020 HC)
hypertension <- c('4079-0.0', '4080-0.0', 'Hypertension')
renal <- c('30670-0.0', '30700-0.0')

# simplify the code
tiss <- 'blood'
arg <- 'severity'
imputed <- TRUE
# note: can repeat the same code, but with exclude_HbA1c = FALSE and vars only containing DM and hba1c
# move things that need manual setting in the same place
# exclude_HbA1c <- TRUE 

# multiple patterns can be applied in one line. (29-Oct-2020 Kai)
# https://community.rstudio.com/t/understanding-the-use-of-str-replace-all-with-multiple-patterns/1849
# dir <- str_replace_all(string = paste0(tiss, '_', Sys.time()), pattern=" |:|-", repl="_") 
dir <- paste0(tiss, '_', format(Sys.time(), "%Y%b%d"))

# make it automatically consider the tissue in study (29-Oct-2020 Kai)
expr_data <- fread(paste0('UKBB_',  tiss, '_predicted_expression_with_pca.csv'), header = T, stringsAsFactors = F)
clinical <- fread(paste0('UKBB_clinical_variables_', ifelse(imputed, 'with', 'without'), '_imputation.csv'))

for(i in c(-1, -3, -7, -8, -9)){
    clinical[clinical == i] <- NA
}

# ---------------------------------------------------------------
# advise to run all variables again (ie change vars) 
# here we have used all the variables selected as covariates, except PCA
# ---------------------------------------------------------------
# new version of selected covariates
traits <- read.csv("selected_covariates_28Oct2020.csv", header = T, stringsAsFactors = F)[[1]]
# ind_PCA = grep('22009', vars) 
vars <- traits[!startsWith(traits, '22009') & (traits != '189-0.0')]
traits <- traits[!startsWith(traits, '22009') & !(traits %in% c('2443-0.1', '2443-1.1'))]
# vars <- read.table('preliminary_study_vars.csv', header = F, stringsAsFactors = F)[[1]]
# selected <- read.table('selected_traits.csv', header = T, stringsAsFactors = F)[[1]]
# statin <- read.table('UKBB_covid_statin_prescription.csv', heaer = T, stringsAsFactors = F)

data <- inner_join(clinical[, c('eid', traits), with = FALSE], expr_data, by = 'eid')
covariates <- data[, -1]
eid <- data[, 1]

# extract outcome variable as vector
responses <- inner_join(clinical[, c('eid', 'U071', 'Severity')], eid, by = 'eid')
if(arg == 'death'){
    outcome <- as.numeric(responses$U071 == 'case')
}else if(arg == 'severity'){
    outcome <- as.numeric(responses$Severity == 'case')
}else{
    stop("arg should be either 'severity' or 'death'!")
}
# change from manually setting to automatic searching (29-Oct-2020 Kai)
# pos <- which(startsWith(colnames(imputed_data), 'ENSG'))[1]
# genetic_vars <- colnames(imputed_data)[pos:ncol(imputed_data)]
# covariates <- imputed_data[, c(traits, genetic_vars), with = F]
# covariates <- sapply(covariates, as.numeric)

is.tuned = F; is_save <- F # HC amended (?to save all results), change it back to False (29-Oct-2020 Kai)
col_names <- c('simes.pval', 'partial.simes.pval', 'pearson.estimate','pearson.pvalue',
               'kendall.estimate','kendall.pvalue', 'spearman.estimate','spearman.pvalue',
               'fisher.pval', 't.test.a.pval', 't.test.b.pval') 

correlation.test.ret <- data.frame(var.name = character(),
                                   simes.pval = double(),
                                   partial.simes.pval = double(),
                                   pearson.est = double(),
                                   pearson.pval = double(),
                                   kendall.est = double(),
                                   kendall.pval = double(),
                                   spearman.est = double(),
                                   spearman.pval = double(),
                                   stringsAsFactors = FALSE)

overall.correlation.ret <- data.frame(var.name = character(),
                                      simes.pval = double(),
                                      partial.simes.pval = double(),
                                      cor.mut.Y.est = double(),
                                      cor.mut.Y.pval = double(),
                                      permutation.pval = double(),
                                      stringsAsFactors = FALSE)

calibration.ret <- data.frame(var.name = character(),
                              calib.mean.est = double(),
                              calib.mean.pval = double(),
                              calib.diff.est = double(),
                              calib.diff.pval = double(),
                              stringsAsFactors = FALSE)

permutate.test.ret <- data.frame(var.name = character(),
                                 var.pval = double(),
                                 risk.pval = double(),
                                 stringsAsFactors = FALSE)


file_prefix <- paste0('covid19_clinical/', dir, '/', arg, '/')
if (!file.exists(file_prefix)){
    dir.create(file_prefix, recursive = T)
}

cat('-------------------------- argument settings --------------------\n')
cat('tissue: ', tiss, '\n')
cat('outcome: ', arg, '\n')
cat('results are saved in ', file_prefix, '\n')

for(tx in vars){
    if(tx == '2443-0.1'){
        treatment <- as.numeric(covariates[['2443-0.0']])
    }else if(tx == '2443-1.1'){
        treatment <- as.numeric(covariates[['2443-1.0']])
    }else{
        treatment <- as.numeric(covariates[, get(tx)])
    }
    
    # i <- which(tx == vars)
    # treatment <- apply(imputed_data[, vars[1:i], with = F], 1, sum)
    included <- (!is.na(treatment) & (treatment >= 0))
    num_excluded <- sum(!included)
    # cat(num_excluded,' NA were excluded.',"\n")
    
    if(length(unique(treatment[included])) == 2){
        binary <- TRUE
    } else {
        binary <- FALSE
    }
	
    # --------------------------------------------------------------------------------
    # HC comment: please note this line may be incorrect if we have continous variable, 
    # as we will force it to zero eg force BMI to zero. 
    # if binary disease trait, can assume NA subjects are controls 
    # If we have enough observations we may remove those NA rows to improve 
    # data accuracy. (29-Oct-2020 Kai)
    # --------------------------------------------------------------------------------
    # if((sum(included) <= 1000) & binary) {
    #     treatment[!included] <- 0
    #     included[!included] <- T
    #     cat(tx,' NA were treated as 0.',"\n")
    # }

    outcome <- outcome[included]
    treatment <- treatment[included]
	
    if(tx %in% smoking_status){
        # HC comment:
        # for example if tx = never smoking, controllng for previous smoking status or current smoking seems strange, 
        # and they are storngly negatively correlated also, perhaps we just exclude 1548 and 1249 from selected traits 
        # eg tx = white population, ie we want to see how the outcome differs between white and non-white,  
        # controlling for whether you are black or asian seems strange.
        # Kai add:
        # We may also need to consider the negative correlation between different ethnic group. For the safety, we need to treat
        # it as one variable.
        X.covariates <- subset(covariates[included,], select = !(colnames(covariates) %in% smoking_status))
    }else if(tx %in% drink_status){
        X.covariates <- subset(covariates[included,], select = !(colnames(covariates) %in% drink_status ))
    }else if(tx %in% ethnic){
        X.covariates <- subset(covariates[included,], select = !(colnames(covariates) %in% ethnic))
    }else if(tx %in% diabetes_exc_HbA1c){
        X.covariates <- subset(covariates[included,], select = !(colnames(covariates) %in% diabetes_exc_HbA1c)) 
    }else if(tx == '2443-0.1'){
        X.covariates <- subset(covariates[included,], select = !(colnames(covariates) %in% diabetes_inc_HbA1c))
    }else if(tx == '2443-1.1'){
        X.covariates <- subset(covariates[included,], select = !(colnames(covariates) %in% diabetes_inc_HbA1c))
    }else if(tx %in% bmi){
        X.covariates <- subset(covariates[included,], select = !(colnames(covariates) %in% bmi))
    }else if(tx %in% hypertension){
        X.covariates <- subset(covariates[included,], select = !(colnames(covariates) %in% hypertension))
    }else if(tx %in% renal){
        X.covariates <- subset(covariates[included,], select = !(colnames(covariates) %in% renal))
    }else{
        X.covariates <- subset(covariates[included,], select = (colnames(covariates) != tx))
    }

    cf.estimator <- ifelse(is.tuned, cf.tuned, cf)

    # split whole dataset into two parts, and the idea of validation is similar to prediction strength.
    pvalues <- split_half_testing(X.covariates, 
                                  outcome, 
                                  treatment, 
                                  binary = binary, 
                                  is_save = is_save, 
                                  num_trees = 5000,
                                  is_tuned = is.tuned,
                                  file_prefix = file_prefix, 
                                  col_names = col_names, 
                                  seed = 123)

    cat(paste0(rep('-', 50), collapse = ''), fill = TRUE)
    cat('variable name:', tx, fill = T)
    cat(num_excluded,' NA were excluded.',"\n")
    cat('pearson correlation pval in trainset:', pvalues[4], fill = T)
    cat('kedall correlation pval in trainset:', pvalues[6], fill = T)
    cat('spearman correlation pval in trainset:', pvalues[8], fill = T)

    # append results to correspoding dataframes
    current_ret <- do.call('c', list(list(tx), as.list(pvalues[1: 8])))
    correlation.test.ret <- append(correlation.test.ret, current_ret)
    # two_sample_test_ret <- do.call('c', list(list(tx), pvalues[9: 11]))
    # double.dataset.test.ret <- append(double.dataset.test.ret, two_sample_test_ret)
                                                                                     
    # *****************************************************************************************
    # fit causal forest on the whole dataset
    # validate fitting with permutating covariates
    # *****************************************************************************************
    tau.forest <- cf.estimator(X.covariates, outcome, treatment, num_trees = 5000)   # run causal forests by default
    tau.prediction <- predict(tau.forest, newdata = NULL, estimate.variance = TRUE, num.threads = 8)
    tau.var <- var(tau.prediction$predictions)
    tau_stats <- compute_stats(tau.prediction)

    # compute zval, pval and ajusted.p
    no.obs <- nrow(X.covariates)
    simes.pval <- simes.test(tau_stats[, 3])
    partial.simes.pval <- simes.partial(floor(no.obs * 0.05), tau_stats[, 3])

    tau_stat <- cbind(eid[[1]][included], outcome, treatment, tau_stats)
    write.csv(tau_stat, file = paste0(file_prefix, tx, '_tau_stats.csv'), quote = F, row.names = F)

    # extract feature importance and save
    varImp <- variable_importance(tau.forest, max.depth = 4)
    varImp.ret <- data.frame(variable = colnames(X.covariates),  varImp)
    write.csv(varImp.ret, paste0(file_prefix, tx, '_varimp.csv'), quote = F, row.names = F) 

    # -----------------------------------------------
    # HC added: extract test_calibration results.
    # -----------------------------------------------
    # cal_obj = test_calibration(tau.forest)
    # test_cal_res = data.frame(cal_obj[1:2, 1:4])
    # write.csv(test_cal_res, paste0(file_prefix, tx, '_test_cal.csv'), quote = F, row.names = F) 

    # test of the calibration with test_calibration from grf 
    # as reported in Github, the pvalue from the method is not accurate. 
    # move the part of extracting and saving test calibration result out of if condition (29-Oct-2020 Kai)
    calibration.fit <- test_calibration(tau.forest)
    cat('mean.pred.estimate of test_calibration:', calibration.fit[1, 1], '; its pval:', calibration.fit[1, 4], fill = TRUE)
    cat('differential.pred.estimate of test_calibration:', calibration.fit[2, 1], '; its pval:', calibration.fit[2, 4], fill = TRUE)
    test_result <- c(calibration.fit[1, c(1, 4)], calibration.fit[2, c(1, 4)])
    record <- do.call('c', list(list(tx), as.list(test_result)))
    calibration.ret <- append(calibration.ret, record)

    # -----------------------------------------------
    # HC: check with ALex on what criteria they used
    # -----------------------------------------------
    # if(simes.pval <= 0.2){ 
    no_sig <- sum(pvalues[c(4, 6, 8)] <= 0.05)

    if(no_sig >= 2){ 
        cor.overall <- cor.test(treatment, outcome, method = 'pearson', alternative = 'greater', use="na.or.complete")
        
        # permutate estimated tau values to validate HET esimation
        Y.hat <- tau.forest[["Y.hat"]]
        W.hat <- tau.forest[["W.hat"]]
        tau.1 <- outcome - Y.hat
        tau.2 <- (treatment - W.hat) * tau.prediction$predictions

        # compute permutated p.value for the mutation
        permutated.p.val <- permutated.pval(tau.1, tau.2)
        cat('p.val by permutating (Y - Y.hat) or (W - W.hat)*tau for ', tx, ':', permutated.p.val, fill = TRUE)

        # save results
        test_result <- c(simes.pval, partial.simes.pval, cor.overall$estimate, cor.overall$p.value, permutated.p.val)
        record <- do.call('c', list(list(tx), as.list(test_result)))
        overall.correlation.ret <- append(overall.correlation.ret, record)

        # valiate hte using by shuffling covariates, and we validate whether the variance of tau after shuffling is likey to 
        # to be less than that observed. 
        # And also, we also use the tau_risk defined by other authors to validate, but in this case, we expect that the tau risk
        # observed is more likey to be smaller tau risk from permutation.
        tau.risk <- assess.explained.tau.fixed.YW.risk(tau.forest, outcome, Y.hat, treatment, W.hat)
        perm.pvals <- adaptive.permutate.covariates.testing(X.covariates, 
                                                            outcome, 
                                                            Y.hat, 
                                                            treatment, 
                                                            W.hat, 
                                                            tau.risk,
                                                            tau.var,  
                                                            is_tuned = is.tuned,
                                                            is_save = is_save, 
                                                            file_prefix = file_prefix, 
                                                            num_trees = 2000, 
                                                            num.strap = 500)
        perm_pval_record <- do.call('c', list(list(tx), as.list(perm.pvals)))
        cat('pval of variance by permutation covariates:', perm.pvals[1], fill = TRUE)
        cat('pval of tau.risk (fixed YW) by permutation covariates:', perm.pvals[2], fill = TRUE)

        permutate.test.ret <- append(permutate.test.ret, perm_pval_record)  
    }
}

write.csv(correlation.test.ret, paste0(file_prefix, 'correlation.test.result.csv'), quote = F, row.names = F)
write.csv(calibration.ret, paste0(file_prefix, 'calibration.result.csv'), quote = F, row.names = F)
write.csv(overall.correlation.ret, paste0(file_prefix, 'overall.correlation.result.csv'), quote = F, row.names = F)
write.csv(permutate.test.ret, paste0(file_prefix, 'permutate.test.result.csv'), quote = F, row.names = F)

