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
        } else if (typeof(x) == "list") {
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
#' @param W_matrix matrix output from `create_tx_matrix` with rows as patient and columns as treatment variables, rows should be ordered the same way as `covar_mat` and `whole_dataset`.
#' @param trainId
#' @param seed default set as 111.
#' @param is_binary whether treatment vector W should be set to binary.
#' @param is_save whether to save all the split observations in SHC.
#' @param is_tuned whether to allow tree self-tuning.
#' @param n_core number of cores to use in paralelle run.
#' @param output_directory file paths of output files, should be created before HTE run.
#' @param skip_perm option to override permutation requirement for quicker run.
#' @param perm_all option to permute all treatment variables 
#' @param random_rep_seed boolean option to set whether different seeds should be generateed for each repeat, by default seeds will be saved for each run to ensure reproducability
#' @return list of data frames including correlation.test.ret, calibration.ret, double.dataset.test.ret, permutate.testing.ret, observed.tau.risk.var.ret, grf.ate.ret, grf.ape.ret, grf.blp.A0.ret. In this order.
#' @export

run.hte <- function(covar_mat,
                    tx_vector,
                    whole_dataset,
                    project,
                    diffCovarTxTypes = FALSE,
                    W_matrix = NULL,
                    trainId,
                    seed = NULL,
                    is_binary = TRUE,
                    is_save = T,
                    save_split = T,
                    is_tuned = F,
                    n_core = 8,
                    output_directory = NULL,
                    skip_perm = FALSE,
                    perm_all = FALSE,
                    random_rep_seed = TRUE,
                    run_blp = TRUE) {
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

    # saving GRF native analysis
    grf.ate.ret <- data.frame(treatment = character(), 
                          est = double(), 
                          std.err = double())

    grf.ape.ret <- data.frame(treatment = character(),
                              est = double(), 
                              std.err = double())

    grf.blp.A0.ret <- data.frame(treatment = character(),
                                 est = double(),
                                 std.error = double(),
                                 t.value = double(),
                                 p.val = double())

    tc_res <- NULL
    ape_res <- NULL
    blp_res <- NULL
    ate_res <- NULL

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

        if (dim(W_matrix)[2] == 1) { # you will get an error if you try to subset via 
            treatment <- W_matrix[,1]
        } else {
            treatment <- W_matrix[,colnames(W_matrix) == tx] # directly retreive the W vector from function input
        }

        if (is_binary) {
            message("Current treatment proportion: ")
            print(table(treatment))
        }

        if (length(unique(treatment)) == 1) {
            print("Not enough obseravtaions for this treatment, skipping.")
            next
        }

        if (!diffCovarTxTypes) X.covariates <- as.matrix(dplyr::select(covar_mat, -all_of(tx)))

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
                                seed = seed,
                                random_rep_seed = random_rep_seed)
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

        if ((simes.pval <= 0.05 & skip_perm == FALSE) | perm_all) {
            print("Performing permutation.")
            cor.overall <- cor.test(covar_mat[, all_of(tx)], Y, method = "pearson", alternative = "greater", use = "na.or.complete")

            # save the result
            pred.ret <- cbind(whole_dataset$donorId, tau_stats)
            colnames(pred.ret) <- c("donorId", "tau.val", "tau.zval", "tau.pval", "tau.p.adjust")
            write.csv(pred.ret, paste0(output_directory, project, "_tau_", tx, ".csv"), quote = F, row.names = F)
            tau.pred.save <- cbind(whole_dataset$donorId, tau.prediction) # saving raw output from pred.causal_forest() to include error estimates
            write.csv(tau.pred.save, file = paste0(file_prefix, "_tau_pred_", tx, ".csv"), row.names = FALSE)

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

            # Add native GRF analysis output

            # average partial effect, should equal ATE for binary unconfounded treatment
            try(cf_ape <- average_partial_effect(tau.forest))
            # best linear prediction of beta0 only
            try(cf_blp_A0 <- best_linear_projection(tau.forest))
            # average treatment effect
            try(cf_ate <- average_treatment_effect(tau.forest))

            if( !all(c(exists("cf_ape"),exists("cf_blp_A0"),exists("cf_ate"))) ) {
                message("Failed to analyse forest.")
                next
            } else if ( any(c( any(is.na(cf_ape)), any(is.na(cf_blp_A0)), any(is.na(cf_ate)) ) ) ) { # evaluate NAs and missing variable separately
                message("Failed to analyse forest.")
                next
            }

            cf_ape <- c(list(tx), as.list(cf_ape))
            grf.ape.ret <- append(grf.ape.ret, cf_ape)

            cf_blp_A0 <- c(list(tx), as.list(cf_blp_A0))
            grf.blp.A0.ret <- append(grf.blp.A0.ret, cf_blp_A0)

            cf_ate <- c(list(tx), as.list(cf_ate))
            grf.ate.ret <- append(grf.ate.ret, cf_ate)

            # best linear prediction conditioned on coviarates
            # reduce covariates until not all of the covairate p values were NAs
            sel_covar <- head(varImp.ret[order(varImp.ret$varImp, decreasing = TRUE),], n = 500)
            A_mat <- X.covariates[, sel_covar$variable] # selected n*500 matrix

            # You will need to manually handle missing values the way GRF handles the missing values
            # see https://grf-labs.github.io/grf/REFERENCE.html#missing-values
            # however, we are not doing splits here so we will need to impute

            A_mat <- knn.impute(A_mat, k = 10)

            try(cf_blp <- best_linear_projection(tau.forest, A = A_mat))
            if (any(class(cf_blp) == "try-error")) next # must check error, it constantly fails with a colname error
            all_na <- all(is.na(cf_blp[,4]))
            reduce <- 500 * (0.9)
            while (all_na){
                sel_covar <- head(varImp.ret[order(varImp.ret$varImp, decreasing = TRUE),], n = reduce)
                if (dim(sel_covar)[1] == 0) {
                    message("BLP failed to plot because all covariates resulted in NAs.")
                }
                A_mat <- X.covariates[, sel_covar$variable] # selected reduced matrix
                cf_blp <- best_linear_projection(tau.forest, A = A_mat)
                all_na <- all(is.na(cf_blp[,4]))
                reduce <- reduce * 0.9
            }

            if (dim(sel_covar)[1] == 0) {
                message("In the process of eliminating covariates s.t. the model isn't all NAs, you have eliminated all covariates. Nothing to see here, moving on.")
            } else {
                message("At least one of the covriates was significant, proceed with blackwards elimination.")
                # Backward elimination of covariates
                elim_blp <- cf_blp[order(cf_blp[,4], decreasing = TRUE),]
                elim_A_mat <- A_mat[, - which(colnames(A_mat) == attr(elim_blp, 'dimnames')[[1]][1])]
                try(elim_blp <- best_linear_projection(tau.forest, A = elim_A_mat))
                if (any(class(elim_blp) == "try-error")) next # must check error, it constantly fails with a colname error
                # some times it gives you both "matrix" and "array" as class, which you will then get the repeated 'Error in dimnames(x) <- dn : length of 'dimnames' [2] not equal to array extent' error
                all_sig <- all(elim_blp[,4] < 0.05)

                while(!all_sig) {
                    # stopping condition when all covariates had been exhausted
                    if (dim(as.matrix(elim_A_mat))[2] == 1) { # at the end of the elimination, when there is only one column left, it will be coerced into a vector and using dim(matrix) alone will return an error
                        message("All covariates were backwards eliminated.")
                        break # exist loop
                    }
                    elim_blp <- elim_blp[order(elim_blp[,4], decreasing = TRUE),]
                    to_b_eliminated <- attr(elim_blp, 'dimnames')[[1]][1]

                    # save to be eliminated covariates
                    to_b_elim_coeff <- elim_blp[1,]
                    to_b_elim_vector <- elim_A_mat[, which(colnames(elim_A_mat) == to_b_eliminated)]

                    if (to_b_eliminated == "(Intercept)") to_b_eliminated <- attr(elim_blp, 'dimnames')[[1]][2] # if the intercept has the highest p value, then eliminate the next in line covariate
                    # message(paste0("Eliminated: ", to_b_eliminated))
                    elim_A_mat <- elim_A_mat[, - which(colnames(elim_A_mat) == to_b_eliminated)]
                    if (dim(as.matrix(elim_A_mat))[2] == 1) { # when only one coviarate is left in the covariate matrix subset, the best_linear_project function in GRF will report a column name error because the subset will be coerced into a colname-less vector
                        elim_A_mat <- as.matrix(elim_A_mat)
                        colnames(elim_A_mat) <- to_b_eliminated
                    }
                    try(elim_blp <- best_linear_projection(tau.forest, A = elim_A_mat))
                    if (any(class(elim_blp) == "try-error")) {
                        all_sig <- TRUE # exit while loop
                        message("Failed to calculate BLP.")
                    } else {
                        all_sig <- all(elim_blp[,4] < 0.05)
                        if (is.na(all_sig)) {
                            # restore previously eliminated covariate
                            elim_A_mat <- cbind(elim_A_mat, to_b_elim_vector)
                            colnames(elim_A_mat)[dim(elim_A_mat)[2]] <- to_b_eliminated
                            try(elim_blp <- best_linear_projection(tau.forest, A = elim_A_mat))
                            all_sig <- TRUE # break while loop
                            message("All coeff in BLP has become NAs.")
                        }
                    }
                }

                print(elim_blp)
                write.csv(elim_blp, paste0(output_directory, project, "_back_elim_blp_", tx, ".csv"), quote = F, row.names = T)
            }
        } else {
            print("Skipping permutation.")
        }

    } # treatment loop
    return(list(correlation.test.ret, 
                calibration.ret, 
                double.dataset.test.ret, 
                permutate.testing.ret, 
                observed.tau.risk.var.ret, 
                grf.ate.ret, 
                grf.ape.ret, 
                grf.blp.A0.ret,
                elim_blp))
} # end of function
