# Functions for running HTE


extract.freq.mutation <- function(dataset, pos, threshold, is_binary = TRUE) {
    # @dataset mutation dataset
    # @pos the position where mutation variables start
    # @threshold threshold for frequent mutation
    if (is_binary == TRUE) {
        mutation.count <- apply(dataset[, pos:ncol(dataset)], 2, function(x) sum(!x == 0))
    } else {
        mutation.count <- apply(dataset[, pos:ncol(dataset)], 2, function(x) sum(x >= 0.75))
    }

    # find mutations with mutation count greater than the defined threshold
    freq.mutations <- colnames(dataset)[pos:ncol(dataset)][mutation.count >= threshold]
    print(paste0("number of mutation will be studied:", sum(mutation.count >= threshold)))
    return(freq.mutations)
}

# if the generic function existis, i.e., append, we should use the same argument names as exist.
# otherwise, we need to create a new one.
# setGeneric("append") # this declear may be omited since append function already exists.
setMethod("append", signature(x = "data.frame", values = "vector"),
    function(x, values) {
        if (ncol(x) != length(values)) {
            stop("dimension of data.frame does not match with the length of values.")
        } else if(typeof(x) == "list") {
            x[nrow(x) + 1, ] <- values
            x
        } else {
            x[nrow(x) + 1, ] <- as.list(values)
            x
        }
    }
)

#' Main function for initiating the HTE analysis with GRF, SHC and permutation
#' 
#' @param covar_mat n by p matrix (X) with n patients and p covariates. All covariates must be numeric, NAs are allowed as of grf (ver 1.2.0) https://grf-labs.github.io/grf/REFERENCE.html#missing-values
#' @param tx_vector string vector of covariate names that will be taken in turns to be the treatment variable (W). All elements in the `tx_vector` should be a subset of the column names of the `covar_mat` unless `diffCovarTxTypes` is set to TRUE.
#' @param whole_dataset n by (p + 2) matrix with n patients and the p covariates plus `donorId` and `outcome`. `whole_dataset` should be ordered the same way as the `covar_mat`
#' @param project string varable to identify the current patient subgroup, e.g. TCGA cancer types or COVID. For naming output paths only.
#' @param diffCovarTxTypes whether the treatment variabels are among the covaraites.
#' @param txdirct binary vector with the same length as `tx_vector` to indicate whether to take the upper or lower quartile as the treatment group.
#' @param trainId
#' @param seed default set as 111.
#' @param is_binary whether treatment vector W should be set to binary.
#' @param is_save whether to save all the split observations in SHC.
#' @param is_tuned whether to allow tree self-tuning.
#' @param thres quartile threshold to select as treatment group if `is_binary` is set to TRUE.
#' @param n_core number of cores to use in paralelle run.
#' @param output_directory file paths of output files, should be created before HTE run.
#' @param skip_perm option to override permutation requirement for quicker run.

run.hte <- function(covar_mat,
                    tx_vector,
                    whole_dataset,
                    project,
                    diffCovarTxTypes = FALSE,
                    txdirct = NULL,
                    trainId,
                    seed = NULL,
                    is_binary = TRUE,
                    is_save = T,
                    save_split = T,
                    is_tuned = F,
                    thres = 0.75,
                    n_core = 8,
                    output_directory = NULL,
                    skip_perm = FALSE) {
    # @covar_mat: covariates matrix (with treatment assignments as well if each of the covariates are taking turns to be analyzed as treatments). Treatment assignments can be binary or continuous.
    # @tx_vector: a vector of variables that will each be used as treatments
    # @whole_dataset: dataframe with outcome, covariates and treatment assignments
    # @is_tuned is consistent within the scope of function, thus it only needs to be set once here.

    correlation.test.ret <- data.frame(gene = character(),
                                       simes.pval = double(),
                                       partial.simes.pval = double(),
                                       pearson.estimate = double(),
                                       pearson.pvalue = double(),
                                       kendall.estimate = double(),
                                       kendall.pvalue = double(),
                                       spearman.estimate = double(),
                                       spearman.pvalue = double(),
                                       stringsAsFactors = FALSE)

    double.dataset.test.ret <- data.frame(gene = character(),
                                          fisher.pval = double(),
                                          t.test.a.pval = double(),
                                          t.test.b.pval = double(),
                                          stringsAsFactors = FALSE)

    calibration.ret <- data.frame(gene = character(),
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

    permutate.testing.ret <- data.frame(gene = character(),
                                        var.pval = double(),
                                        fixed.YW.risk.pval = double(),
                                        stringsAsFactors = FALSE)

    observed.tau.risk.var.ret <- data.frame(gene = character(),
                                            var = double(),
                                            fixed.YW.risk = double(),
                                            stringsAsFactors = FALSE)



    no.obs <- dim(covar_mat)[1]
    no.obs.train <- length(trainId)
    no.obs.test <- no.obs - no.obs.train

    Y <- whole_dataset$outcome

    cf.estimator <- ifelse(is_tuned, cf.tuned, cf)

    i <- 1
    for (tx in tx_vector) {

        print(paste0(c("#", rep("-", 40), " begin a new treatment ", rep("-", 40)), collapse = ""))
        print(paste0("Processing ", i, " of ", length(tx_vector), " genes."))
        i <- i + 1

        # since treatment variable is 0 or 1, if the feature considered is binary, then no transformation is neeeded;
        # otherwise, we set values of the feature greater than specific quantile, say 0.75, to 1

        treatment <- covar_mat[, tx] # it has been confirmed that grf can deal with continuous treatment variable

        treatment <- assign_tx(binary = is_binary, upperQ = txdirct[[tx]], thres = thres, treatment = treatment) # Use the assign_tx function to convert treatment vector

        # Check if we have enough tx observations
        if (!is_binary) {
            if (length(unique(treatment)) == 1) {
                print("Not enough obseravtaions for this treatment, skipping.")
                next
            }
        }

        if (!diffCovarTxTypes) X.covariates <- as.matrix(dplyr::select(covar_mat, -tx))

        # print(paste0(c('#', rep('-', 40), ' begin a new treatment ', rep('-', 40)), collapse = ''))

        # split whole dataset into two parts, and the idea of validation is similar to prediction strength.

        col_names <- c("simes.pval", "partial.simes.pval", "pearson.estimate",
                        "pearson.pvalue", "kendall.estimate", "kendall.pvalue",
                        "spearman.estimate", "spearman.pvalue", "fisher.pval",
                        "t.test.a.pval", "t.test.b.pval")

        #test
        file_prefix <- paste0(output_directory, project, "_", tx)

        print(file_prefix)

        print("Performing split half.")
        pvalues <- try(split_half_testing(X.covariates, Y,
                                treatment,
                                binary = is_binary,
                                is_save = is_save,
                                save_split = save_split,
                                varimp_names = colnames(X.covariates),
                                is_tuned = is_tuned,
                                file_prefix = file_prefix,
                                col_names = col_names,
                                seed = seed)
            )

        if (class(pvalues) == "try-error") next

        cat("Treatment name: ", tx, fill = T)
        cat("Fisher extact test pval in trainset: ", pvalues[9], fill = T)
        cat("pearson correlation pval in trainset: ", pvalues[4], fill = T)
        cat("kedall correlation pval in trainset: ", pvalues[6], fill = T)
        cat("spearman correlation pval in trainset: ", pvalues[8], fill = T)

        # append results to correspoding dataframes
        current_ret <- do.call("c", list(list(tx), as.list(pvalues[1: 8])))
        correlation.test.ret <- append(correlation.test.ret, current_ret)
        two_sample_test_ret <- do.call("c", list(list(tx), pvalues[9: 11]))
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

        print(paste0("simes.pval is ", simes.pval))

        if (simes.pval <= 0.05 & skip_perm == FALSE) {
            print("Performing permutation.")
            cor.overall <- cor.test(covar_mat[, tx], Y, method = "pearson", alternative = "greater", use = "na.or.complete")

            # save the result
            pred.ret <- cbind(whole_dataset$donorId, tau_stats)
            colnames(pred.ret) <- c("donorId", "tau.val", "tau.zval", "tau.pval", "tau.p.adjust")
            write.csv(pred.ret, paste0(output_directory, project, "_tau_", tx, ".csv"), quote = F, row.names = F)

            # permutate estimated tau values to validate HET esimation
            Y.hat <- tau.forest[["Y.hat"]]
            W.hat <- tau.forest[["W.hat"]]
            tau.1 <- Y - Y.hat
            tau.2 <- (treatment - W.hat) * tau.prediction$predictions

            # compute permutated p.value for the mutation
            permutated.p.val <- permutated.pval(tau.1, tau.2)
            print(paste0("p.val by permutating (Y - Y.hat) or (W - W.hat)*tau for ", tx, ":", permutated.p.val))

            # test of the calibration with test_calibration from grf
            # as reported in Github, the pvalue from the method is not accurate.
            calibration.fit <- test_calibration(tau.forest)
            print(paste0("mean.pred.estimate of test_calibration:", calibration.fit[1, 1], "; its pval:", calibration.fit[1, 4]))
            print(paste0("differential.pred.estimate of test_calibration:", calibration.fit[2, 1], "; its pval:", calibration.fit[2, 4]))

            # save results
            test_result <- c(calibration.fit[1, c(1, 4)], calibration.fit[2, c(1, 4)])
            test_result <- c(simes.pval, partial.simes.pval, cor.overall$estimate, cor.overall$p.value, permutated.p.val, test_result)
            record <- do.call("c", list(list(tx), as.list(test_result)))
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
                                                                is_tuned = is_tuned,
                                                                is_save = is_save,
                                                                file_prefix = file_prefix,
                                                                num_trees = 1000,
                                                                num.strap = 500)

            perm_pval_record <- do.call("c", list(list(tx), as.list(perm.pvals)))
            perm_var_risk_record <- do.call("c", list(list(tx), as.list(c(tau.var, fixed.YW.tau.risk))))

            permutate.testing.ret <- append(permutate.testing.ret, perm_pval_record)
            observed.tau.risk.var.ret <- append(observed.tau.risk.var.ret, perm_var_risk_record)

            print(paste0("pval of variance by permutation covariates:", perm.pvals[1]))
            print(paste0("pval of tau.risk (fixed YW) by permutation covariates:", perm.pvals[2]))

            # extract feature importance and save
            varImp <- variable_importance(tau.forest,  max.depth = 4)
            varImp.ret <- data.frame(variable = colnames(X.covariates),  varImp)
            write.csv(varImp.ret, paste0(output_directory, project, "_varimp_", tx, ".csv"), quote = F, row.names = F)
        } else {
            print("Skipping permutation.")
        }


    }
    return(list(correlation.test.ret, calibration.ret, double.dataset.test.ret, permutate.testing.ret, observed.tau.risk.var.ret))
}
