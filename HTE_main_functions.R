# Functions for running HTE


extract.freq.mutation <- function(dataset, pos, threshold, is.binary = TRUE){
    # @dataset mutation dataset
    # @pos the position where mutation variables start
    # @threshold threshold for frequent mutation
    if(is.binary == TRUE){
        mutation.count <- apply(dataset[, pos:ncol(dataset)], 2, function(x) sum(!x == 0))
    }else{
        mutation.count <- apply(dataset[, pos:ncol(dataset)], 2, function(x) sum(x >= 0.75))
    }

    # find mutations with mutation count greater than the defined threshold
    freq.mutations <- colnames(dataset)[pos:ncol(dataset)][mutation.count >= threshold]
    print(paste0('number of mutation will be studied:', sum(mutation.count >= threshold)))
    return(freq.mutations)
}

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

run.hte <- function(covar_mat, tx_vector, whole_dataset, project, covar_type = NULL, trainId, seed = NULL, is.binary = TRUE, is_save = T, save_split = T, is.tuned = F, thres = 0.75, n_core = 8, output_directory = NULL){
    # @covar_mat: covariates matrix (with treatment assignments as well if each of the covariates are taking turns to be analyzed as treatments). Treatment assignments can be binary or continuous.
    # @tx_vector: a vector of variables that will each be used as treatments
    # @whole_dataset: dataframe with outcome, covariates and treatment assignments
    # @is.tuned is consistent within the scope of function, thus it only needs to be set once here.
    
    correlation.test.ret <- data.frame(mutName = character(),
                                       simes.pval = double(),
                                       partial.simes.pval = double(),
                                       pearson.estimate = double(),
                                       pearson.pvalue = double(),
                                       kendall.estimate = double(),
                                       kendall.pvalue = double(), 
                                       spearman.estimate = double(),
                                       spearman.pvalue = double(),
                                       stringsAsFactors = FALSE)

    double.dataset.test.ret <- data.frame(mutName = character(),
                                          fisher.pval = double(),
                                          t.test.a.pval = double(),
                                          t.test.b.pval = double(),
                                          stringsAsFactors = FALSE)

    calibration.ret <- data.frame(mutName = character(),
                                  simes.pval = double(),
                                  simes.partial = double(),
                                  cor.mut.Y.estimate = double(),
                                  cor.mut.Y.pvalue = double(),
                                  permutation.pval = double(),
                                  mean.pred.estimate = double(),
                                  mean.pred.pval = double(),
                                  diff.pred.estimate = double(),
                                  diff.pred.pval = double(),
                                  stringsAsFactors = FALSE)

    permutate.testing.ret <- data.frame(mutName = character(),                                       
                                        var.pval = double(),
                                        fixed.YW.risk.pval = double(),
                                        stringsAsFactors = FALSE)

    observed.tau.risk.var.ret <- data.frame(mutName = character(),                                       
                                            var = double(),
                                            fixed.YW.risk = double(),
                                            stringsAsFactors = FALSE)
    
    
    
    no.obs <- dim(covar_mat)[1]
    no.obs.train <- length(trainId)
    no.obs.test <- no.obs - no.obs.train

    Y <- whole_dataset$outcome

    cf.estimator <- ifelse(is.tuned, cf.tuned, cf)
    
    i = 1
    for(tx in tx_vector){
        
        # since treatment variable is 0 or 1, if the feature considered is binary, then no transformation is neeeded;
        # otherwise, we set values of the feature greater than specific quantile, say 0.75, to 1

        treatment <- covar_mat[, tx] # it has been confirmed that grf can deal with continuous treatment variable


        X.covariates <- as.matrix(dplyr::select(covar_mat, -tx))
        
        print(paste0(c('#', rep('-', 40), ' begin a new treatment ', rep('-', 40)), collapse = ''))
        
        print(paste0("Processing ", i, " of ", length(tx_vector), " genes."))
        i = i+1
        
        if(is.binary){
            if (covar_type == "mutation") {
                treatment <- as.numeric(treatment != 0) # only for mutation
                if (length(unique(treatment)) == 1) {
                    print("Gene mutation distribution too sparse, skipping.")
                    next
                }
                
            } else if (covar_type == "expression") {
               
               # Read corresponding FC in tumor
                    DEGs = read.csv(paste0("./tables/", project, "_DEGtable.csv"))
               # determine if over expressed or under expressed in tumor
                    DEmeters = dplyr::filter(DEGs, X==tx)
               # set threshold based on gene behavior in tumors
                    if (DEmeters$logFC > 0){
                        # we take UQ
                        treatment = as.numeric(treatment > quantile(treatment, 0.75))
                    } else {
                        # we take LQ
                        treatment = as.numeric(treatment < quantile(treatment, 0.25))
                        }
            
                    
                    # Safeguarding against uniform treatment assignment
                    if (length(unique(treatment)) == 1 | sum(treatment != 0) < length(treatment)*0.1) {
                        print("Gene expression distribution too sparse, skipping.")
                        next
                    }
                    
            } else if (covar_type=="UQ"){
                    treatment = as.numeric(treatment > quantile(treatment, 0.75))
                    if (length(unique(treatment)) == 1 | sum(treatment) < length(treatment)*0.1) {
                        print("Gene expression distribution too sparse, skipping.")
                        next
                }
            }
                
                
                 else {
                    treatment <- as.numeric(treatment)
                    }
        }
        # split whole dataset into two parts, and the idea of validation is similar to prediction strength.
        
        col_names <- c('simes.pval', 'partial.simes.pval', 'pearson.estimate','pearson.pvalue',
               'kendall.estimate','kendall.pvalue', 'spearman.estimate','spearman.pvalue',
               'fisher.pval', 't.test.a.pval', 't.test.b.pval')

        #test
        file_prefix = paste0(output_directory, project, "_", tx)

        print(file_prefix)
        
        print("Performing split half.")
        pvalues <- split_half_testing(X.covariates, Y, 
                                      treatment, 
                                      binary = is.binary, 
                                      is_save = is_save, 
                                      save_split = save_split, 
                                      is_tuned = is.tuned,
                                      file_prefix = file_prefix,
                                      col_names = col_names, 
                                      seed = seed) 


        cat('Treatment name:', tx, fill = T)
        cat('Fisher extact test pval in trainset:', pvalues[9], fill = T)
        cat('pearson correlation pval in trainset:', pvalues[4], fill = T)
        cat('kedall correlation pval in trainset:', pvalues[6], fill = T)
        cat('spearman correlation pval in trainset:', pvalues[8], fill = T)
        
        # append results to correspoding dataframes
        current_ret <- do.call('c', list(list(tx), as.list(pvalues[1: 8])))
        correlation.test.ret <- append(correlation.test.ret, current_ret)
        two_sample_test_ret <- do.call('c', list(list(tx), pvalues[9: 11]))
        double.dataset.test.ret <- append(double.dataset.test.ret, two_sample_test_ret)
                                                                                         
        # *****************************************************************************************
        # fit causal forest on the whole dataset
        # validate fitting with permutating covariates
        # *****************************************************************************************
        print("Fitting CF on the whole dataset.")        
        tau.forest <- cf.estimator(X.covariates, Y, treatment)   # run causal forests by default
        tau.prediction <- predict(tau.forest, newdata = NULL, estimate.variance = TRUE, num.threads = n_core)
        tau.var <- var(tau.prediction$predictions)

        # compute zval, pval and ajusted.p
        tau_stats <- compute_stats(tau.prediction)
        simes.pval <- simes.test(tau_stats[, 3])
        partial.simes.pval <- simes.partial(floor(no.obs * 0.05), tau_stats[, 3])

        if(simes.pval <= 0.05) { 
            print("Performing permutation.")
            cor.overall <- cor.test(covar_mat[, tx], Y, method = 'pearson', alternative = 'greater', use="na.or.complete")
            
            # save the result
            pred.ret <- cbind(whole_dataset$donorId, tau_stats) 
            colnames(pred.ret) <- c('donorId', 'tau.val', 'tau.zval', 'tau.pval', 'tau.p.adjust')
            write.csv(pred.ret, paste0(output_directory, project, '_tau_', tx, '.csv'), quote = F, row.names = F)

            # permutate estimated tau values to validate HET esimation
            Y.hat <- tau.forest[["Y.hat"]]
            W.hat <- tau.forest[["W.hat"]]
            tau.1 <- Y - Y.hat
            tau.2 <- (treatment - W.hat) * tau.prediction$predictions

            # compute permutated p.value for the mutation
            permutated.p.val <- permutated.pval(tau.1, tau.2)
            print(paste0('p.val by permutating (Y - Y.hat) or (W - W.hat)*tau for ', tx, ':', permutated.p.val))

            # test of the calibration with test_calibration from grf 
            # as reported in Github, the pvalue from the method is not accurate. 
            calibration.fit <- test_calibration(tau.forest)
            print(paste0('mean.pred.estimate of test_calibration:', calibration.fit[1, 1], '; its pval:', calibration.fit[1, 4]))
            print(paste0('differential.pred.estimate of test_calibration:', calibration.fit[2, 1], '; its pval:', calibration.fit[2, 4]))

            # save results
            test_result <- c(calibration.fit[1, c(1, 4)], calibration.fit[2, c(1, 4)])
            test_result <- c(simes.pval, partial.simes.pval, cor.overall$estimate, cor.overall$p.value, permutated.p.val, test_result)
            record <- do.call('c', list(list(tx), as.list(test_result)))
            calibration.ret <- append(calibration.ret, record)

            # valiate hte using by shuffling covariates, and we validate whether the variance of tau after shuffling is likey to 
            # to be less than that observed. 
            # And also, we also use the tau_risk defined by other authors to validate, but in this case, we expect that the tau risk
            # observed is more likey to be smaller tau risk from permutation.
            fixed.YW.tau.risk <- assess.explained.tau.fixed.YW.risk(tau.forest, Y, Y.hat, treatment, W.hat)

            perm.pvals <- adaptive.permutate.covariates.testing(X.covariates, 
                                                                Y, Y.hat, 
                                                                treatment, 
                                                                W.hat, 
                                                                fixed.YW.tau.risk,
                                                                tau.var,  
                                                                is_tuned = is.tuned,
                                                                is_save = is_save, 
                                                                file_prefix = file_prefix, 
                                                                num_trees = 1000, 
                                                                num.strap = 500)
            
            perm_pval_record <- do.call('c', list(list(tx), as.list(perm.pvals)))
            perm_var_risk_record <- do.call('c', list(list(tx), as.list(c(tau.var, fixed.YW.tau.risk))))

            permutate.testing.ret <- append(permutate.testing.ret, perm_pval_record)    
            observed.tau.risk.var.ret <- append(observed.tau.risk.var.ret, perm_var_risk_record)    

            print(paste0('pval of variance by permutation covariates:', perm.pvals[1]))
            print(paste0('pval of tau.risk (fixed YW) by permutation covariates:', perm.pvals[2]))

            # extract feature importance and save
            varImp <- variable_importance(tau.forest,  max.depth = 4)
            varImp.ret <- data.frame(variable = colnames(X.covariates),  varImp)
            write.csv(varImp.ret, paste0(output_directory, project, '_varimp_', tx, '.csv'), quote = F, row.names = F) 
        } else {
            print("Skipping permutation.")
        }
        
        
    }
    return(list(correlation.test.ret, calibration.ret, double.dataset.test.ret, permutate.testing.ret, observed.tau.risk.var.ret))
}
