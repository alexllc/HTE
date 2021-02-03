# Script for running testing whether treatment group proportion affects results significantly

## Run from base directory
setwd("~/project/HTE/")

## source bin scripts
bin_ls = list.files("./bin")
for (bin in bin_ls){
    source(paste0("./bin/", bin))
}

## Set parameters for this run
cancer_type = "BRCA"
endpt = "OS"
output_directory = "./exp/2020-12-11_tx_gruop_size_test/tmp/interactive_out/"
result_directory = "./exp/2020-12-11_tx_gruop_size_test/res/"
is_tuned = FALSE
is_save = FALSE
n_core = 40
seed = 111

## Prepare clinical dataframe
clinical = fetch_clinical_data(cancer_type, outParam = endpt, imputeMethod = "simple", outUnitDays2Month = TRUE, discard = c("type", "tumor_status"))
clinical = mk_id_rownames(clinical)

## Prepare expression dataframe
exp = fetch_exp_data(cancer_type, scale = FALSE, primaryTumorOnly = TRUE, formatPatient = TRUE)

## Subset patients and settings for HTE
# Check the list of common patients across three data frames
common_pat = rownames(clinical)[rownames(clinical) %in% rownames(exp)]
message(paste0("Number of patients with expression, mutation and clinical entries is: ", length(common_pat)))

# Outcome vector for causal forest
Y = clinical[common_pat, "outcome"]
X = exp[common_pat,]
wds <- merge(X, clinical[common_pat,], by = 0, suffix = c("",""))
colnames(wds)[colnames(wds) == "Row.names"]  <- "donorId"
X <- dplyr::select(wds, -c(outcome, donorId))

# Variables copied from HTE_main_functions.R
obsNumber <- dim(X)[1]
trainId <- sample(1: obsNumber, floor(obsNumber/2), replace = FALSE)
registerDoParallel(10)
no.obs <- dim(X)[1]
no.obs.train <- length(trainId)
no.obs.test <- no.obs - no.obs.train

# empty dataframes for initiating results
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

calibration_summary <- data.frame(  txProp = double(),
                                    simesPval = double(),
                                    meanPredEst = double(),
                                    meanPredPval = double(),
                                    diffPredEst = double(),
                                    diffPredPval = double())

permutate_summary <- data.frame(txProp = double(),
                                varPval = double(),
                                fixYWRiskPval = double())

# Import DEA results to determine treatment directions
DEA_tbl <- read.csv(paste0("./dat/tables/", cancer_type, "_DEGtable.csv"))
rownames(DEA_tbl) <- DEA_tbl$X

# Will no longer test only significant genes ran in the past, but all genes that are present in the DEA analysis result
sig_ls <- colnames(exp)[colnames(exp) %in% DEA_tbl$X]

prop_ls <- c(0.01, 0.02, 0.05, 0.1, 0.25, 0.5)
prop_names <- c("one", "two", "five", "ten", "quarter", "half")

for(prop in prop_ls) {
    message(rep("=", 80))
    message(paste0("Proportion of patients used :", prop))
    message(rep("=", 80))
    for(tx in sig_ls) {
        message(paste0("Processing gene: ", tx, " ( ", which(sig_ls %in% tx), " out of ", length(sig_ls), " )"))
        if (DEA_tbl[tx,"logFC"] > 0) {
            treatment <- as.numeric(X[,colnames(X)==tx] > quantile(X[,colnames(X)==tx], 1 - prop)) # in the actual result we chose, we will use the opposite direction rather than the direction along with DEA result
        } else {
           treatment <- as.numeric(X[,colnames(X)==tx] < quantile(X[,colnames(X)==tx], prop))
        }
        
        file_prefix <- paste0(output_directory, cancer_type, "_", tx)

        X.covariates <- dplyr::select(X, -c(tx))
        # The following code is a direct copy from the HTE_main_functions, run.hte, I don't want make a duplicate function of this part.

        cf.estimator <- ifelse(is_tuned, cf.tuned, cf)



        tau.forest <- cf.estimator(X.covariates, Y, treatment)   # run causal forests by default
        tau.prediction <- predict(tau.forest, newdata = NULL, estimate.variance = TRUE, num.threads = n_core)
        tau.var <- var(tau.prediction$predictions)

        # compute zval, pval and ajusted.p
        tau_stats <- compute_stats(tau.prediction)
        simes.pval <- simes.test(tau_stats[, 3])
        if (is.na(simes.pval)) {
            message("Invalid Simes p values.")
            next
        } else if (simes.pval > 0.05) {
           message("Non-significant Simes p value.")
           next
        }
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

    assign(paste0(prop_names[which(prop_ls %in% prop)], "_calibrate"), calibration.ret)
    assign(paste0(prop_names[which(prop_ls %in% prop)], "_permutate"), permutate.testing.ret)
    assign(paste0(prop_names[which(prop_ls %in% prop)], "_obs_tau_risk"), observed.tau.risk.var.ret)

    write.csv(get(paste0(prop_names[which(prop_ls %in% prop)], "_calibrate")), 
            file = paste0(output_directory, prop, "_calibrate.csv"), row.names = FALSE)
    write.csv(get(paste0(prop_names[which(prop_ls %in% prop)], "_permutate")), 
            file = paste0(output_directory, prop, "_permutate.csv"), row.names = FALSE)
    write.csv(get(paste0(prop_names[which(prop_ls %in% prop)], "_obs_tau_risk")), 
            file = paste0(output_directory, prop, "_obs_tau_risk.csv"), row.names = FALSE)

    calibration_rc <- do.call("c", list(list(prop), as.list(c(simes.test(calibration.ret[,2]), mean(calibration.ret[,7]), simes.test(calibration.ret[,8]), mean(calibration.ret[,9]), simes.test(calibration.ret[,10])) )))
    permutate_rc <- do.call("c", list(list(prop), as.list( c(simes.test(permutate.testing.ret[,2]), simes.test(permutate.testing.ret[,3])) )))

    calibration_summary <- append(calibration_summary, calibration_rc)
    permutate_summary <- append(permutate_summary, permutate_rc)
}

write.csv(calibration_summary, file = paste0(result_directory, cancer_type, "_no_select_calibration_simes_pval_summary.csv"), row.names = FALSE)
write.csv(permutate_summary, file = paste0(result_directory, cancer_type, "_no_select_permutate_simes_pval_summary.csv"), row.names = FALSE)