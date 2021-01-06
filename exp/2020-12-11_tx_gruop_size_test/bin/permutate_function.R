

        tau.forest <- cf.estimator(X.covariates, Y, treatment)   # run causal forests by default
        tau.prediction <- predict(tau.forest, newdata = NULL, estimate.variance = TRUE, num.threads = n_core)
        tau.var <- var(tau.prediction$predictions)

        # compute zval, pval and ajusted.p
        tau_stats <- compute_stats(tau.prediction)
        simes.pval <- simes.test(tau_stats[, 3])
        partial.simes.pval <- simes.partial(floor(no.obs * 0.05), tau_stats[, 3])

        print(paste0("simes.pval is ", simes.pval))

        print("Performing permutation.")
        cor.overall <- cor.test(X[, colnames(X) == tx], Y, method = "pearson", alternative = "greater", use = "na.or.complete")

        # save the result
        pred.ret <- cbind(wds$donorId, tau_stats)
        colnames(pred.ret) <- c("donorId", "tau.val", "tau.zval", "tau.pval", "tau.p.adjust")
        write.csv(pred.ret, paste0(output_directory, cancer_type, "_tau_", tx, ".csv"), quote = F, row.names = F)

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
    }
    assign(paste0(prop, "_calibrate"), calibration.ret)
    assign(paste0(prop, "_permutate"), permutate.testing.ret)
    assign(paste0(prop, "_obs_tau_risk"), observed.tau.risk.var.ret)

    calibration_rc <- do.call("c", list(list(prop), as.list(c(simes.test(calibration.ret[,2]), mean(calibration.ret[,7]), simes.test(calibration.ret[,8]), mean(calibration.ret[,9]), simes.test(calibration.ret[,10])) )))
    permutate_rc <- do.call("c", list(list(prop), as.list( c(simes.test(permutate.testing.ret[,2]), simes.test(permutate.testing.ret[,3])) )))

    calibration_summary <- append(calibration_summary, calibration_rc)
    permutate_summary <- append(permutate_summary, permutate_rc)
