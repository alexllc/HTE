# some methods for validate causal models
# for bootstrap testing refer to "a comparison of methods for model selection when estimating individual treatment effects"
# permutate.covariates.testing is to test whether covariates contribute to variation in heterogeneous treatment effect.
# for permutated.pval ref. equation 9 to the paper above, we permute y_i - m_hat(x_i) to check our fitting

library(DescTools)


early.stop <- function(pvals){
    is_stop <- (sum(pvals > 0.05) == length(pvals))
    is_stop
}

get.lower.BinomCI <- function(perm.val, observed.val, min_perm, conf.level, BinomCI_method){
    CI <- BinomCI(x = sum(perm.val > observed.val), 
                  n = min_perm, 
                  conf.level = conf.level, 
                  method = BinomCI_method)
    CI[2]
}

extract_binom_pval <- function(est){
    binom_obj <- binom.test(x = sum(est > 0), 
                            n = 20 , 
                            p = 0.5,
                            alternative = "greater",
                            conf.level = 0.95)
    return(binom_obj$p.value)
}

simes.test <- function(x, returnstat = FALSE){
    r = rank(x,  ties.method = "random")
    t = min(length(x) * x / r)
    if (returnstat) c(t, t) else t
}

simes.partial <- function(u, pvec)  {  # equation 2 in the above reference
    n = length(pvec)
    pvec = pvec[order(pvec)]  # sort pval in ascending order

    corr.p = numeric(n-u+1)

    for (i in 1: (n - u + 1) ){
        corr.p[i] = (n-u+1)/i  * pvec[u-1+i]
    }

    partial.conj.p = min(corr.p)
    return(partial.conj.p)
}

Simes.partial.allsubset <- function(u, pvec)  {  ##equation 2 in the above reference
    n <- length(pvec)
    # pvec = pvec[order(pvec)]  ##sort pval in ascending order
  
    setsize <- n - u + 1
    subset_mat <- subsets(n = n, r = setsize, pvec)
    all_simes_res <- apply(subset_mat, 1, simes.test) 
    return(list(all_simes_res = all_simes_res, Simes_allsubset = max(all_simes_res)))
}

harmonMeanP.partial.allsubset <- function(u, pvec)  {  
  n <- length(pvec)
  #pvec = pvec[order(pvec)]  ##sort pval in ascending order
  
  setsize <- n - u + 1
  subset_mat = subsets(n=n, r=setsize, pvec)
  
  all_res = apply(subset_mat, 1, p.hmp) 
  return(list(all_res=all_res, harmonMeanP_allsubset = max(all_res) )  )

}

Fisher.partial <- function(u, pvec) {
   n = length(pvec)
   pvec = pvec[order(pvec)]  ##sort pval in ascending order
   chi_sq_thres = (-2)*sum(log(pvec[u:n] )) 
   Fisher_partial_conjunction_p <- pchisq(q =  chi_sq_thres, 
                                          df = 2 * (n - u + 1) , 
                                          ncp = 0, 
                                          lower.tail = FALSE, 
                                          log.p = FALSE)
            
   return(Fisher_partial_conjunction_p)
}

stratified.permutation <- function(cluster_id){
    cluters <- levels(cluster_id)
    obs <- seq(length(cluster_id))

    samp <- rep(0, length(cluster_id))
    for(cls in cluters){
        cls_idx <- obs[cluster_id == cls]
        samp[cls_idx] <- sample(cls_idx, length(cls_idx), replace = F)
    }
    return(samp)
}

assess.explained.tau.fixed.YW.risk <- function(fitted.obj, Y, Y.hat, treatment, W.hat, constant.tau = NULL){
    # calculate tau-risk_R
    tau.pred <- predict(fitted.obj, newdata = NULL, estimate.variance = TRUE, num.threads = 8)
    tau.1 <- Y - Y.hat

    if (is.null(constant.tau)){
        mean.tau <- (treatment - W.hat) * mean(tau.pred$predictions)
        tau.2 <- (treatment - W.hat) * tau.pred$predictions
    }else {
        mean.tau <- (treatment - W.hat) * mean(constant.tau)
        tau.2 <- (treatment - W.hat) * constant.tau
    }
    
    risk.total <- (tau.1 - mean.tau)^2
    risk.unexplained <- (tau.1 - tau.2)^2
    risk.explained <- mean(risk.total - risk.unexplained)
    return(risk.explained)
}

split_half_testing <- function(covariates, Y, 
                               treatment, 
                               binary = F, 
                               is_save = T, 
                               save_split = T, 
                               file_prefix = NULL, 
                               col_names = NULL, 
                               is_tuned = F, 
                               no_repeats = 10, 
                               seed = NULL, # seed will be set within run.hte call 
                               n_core = 8){

    # my added line 20 Sep 2019
    cf.estimator <- ifelse(is_tuned, cf.tuned, cf)

    # note: if is_save is TRUE, then file_prefix and col_names must be provided.
    # when is_save turns on, results for every repeat will be saved. 
    # if save_split is true, then observation level result will be stored.
    # Here file names is a file prefix for saving
    if((is_save | save_split) && is.null(file_prefix)) stop('if is_save or save_split is T, then file_prefix must be provided.')
    
    no.obs <- dim(covariates)[1]
    #correlation_matrix <- foreach(i = seq(no_repeats), .combine = 'rbind', .options.multicore = list(set.seed = T)) %dopar%  {

    # For loop
    correlation_matrix = NULL
    for (i in seq(no_repeats)) {
        observation_result.a <- matrix(0, nrow = no.obs, ncol = 4)
        observation_result.b <- matrix(0, nrow = no.obs, ncol = 4)

        if(binary){
            treat <- seq(no.obs)[treatment == 1]
            control <- seq(no.obs)[treatment == 0]
            sampled_treat <- floor(length(treat)/2)
            trainId <- c(sample(treat, sampled_treat), sample(control, floor(no.obs/2) - sampled_treat))
        }else{
            trainId <- sample(1: no.obs, floor(no.obs/2), replace = FALSE)
        }

        no.obs.train <- length(trainId)
        no.obs.test <- no.obs - no.obs.train

        tau.forest.train <- cf.estimator(covariates[trainId, ], Y[trainId], treatment[trainId], seed = seed)   
        tau.pred.train.a <- predict(tau.forest.train, estimate.variance = T, num.threads = n_core)
        tau.pred.test.a <- predict(tau.forest.train, newdata = covariates[-trainId, ], estimate.variance = T, num.threads = n_core)

        tau.forest.test <- cf.estimator(covariates[-trainId, ], Y[-trainId], treatment[-trainId], seed = seed) 
        tau.pred.train.b <- predict(tau.forest.test, estimate.variance = T, num.threads = n_core)
        tau.pred.test.b <- predict(tau.forest.test, newdata = covariates[trainId, ], estimate.variance = T, num.threads = n_core)

        # compute z-score, pvalues, and ajusted.pvalues
        tau.a.train.stats <- compute_stats(tau.pred.train.a)
        tau.a.stats <- compute_stats(tau.pred.test.a)
        tau.b.train.stats <- compute_stats(tau.pred.train.b)
        tau.b.stats <- compute_stats(tau.pred.test.b)

        if(save_split){
            observation_result.a[trainId,] <- tau.a.train.stats
            observation_result.a[-trainId,] <- tau.a.stats
            write.csv(observation_result.a, file = paste0(file_prefix, project,  '_observation_', i, '_result_a.csv'))
            observation_result.b[-trainId,] <- tau.b.train.stats 
            observation_result.b[trainId,] <- tau.b.stats
            write.csv(observation_result.b, file = paste0(file_prefix, project, '_observation_', i,'_result_b.csv'))
        }

        simes_pval.a <- simes.test(tau.a.stats[, 3])
        simes_pval.b <- simes.test(tau.b.stats[, 3])

        partial_simes_pval.a <- simes.partial(floor(no.obs.test * 0.05), tau.a.stats[, 3])
        partial_simes_pval.b <- simes.partial(floor(no.obs.train * 0.05), tau.b.stats[, 3])

        # check the correlation between two predictions from two datasets
        test_res.a <- correlation_test(tau.a.train.stats[, 1], tau.b.stats[, 1], methods = c('pearson', 'kendall', 'spearman'))
        test_res.b <- correlation_test(tau.a.stats[, 1], tau.b.train.stats[, 1], methods = c('pearson', 'kendall', 'spearman'))

        fisher.pval.a <- fisher.exact.test(tau.a.train.stats[, 3], tau.b.stats[, 3])
        fisher.pval.b <- fisher.exact.test(tau.a.stats[, 3], tau.b.train.stats[, 3])

        t.test.pval.a <- quantile.t.test(tau.a.train.stats[, 1], tau.b.stats[, 1])
        t.test.pval.b <- quantile.t.test(tau.a.stats[, 1], tau.b.train.stats[, 1]) 
        
        correlation_rslt <- rbind(c(simes_pval.a, partial_simes_pval.a, test_res.a, fisher.pval.a, t.test.pval.a), 
                                  c(simes_pval.b, partial_simes_pval.b, test_res.b, fisher.pval.b, t.test.pval.b))
        correlation_matrix = rbind(correlation_matrix, correlation_rslt)
    }

    if(is_save){
        colnames(correlation_matrix) <- col_names 
        write.csv(correlation_matrix, file = paste0(file_prefix, project, '_split_half.csv'), row.names = F, quote = F)
    }   
    
    # change to partial_simes_pval 20190929
    aggregated_rslt <- sapply(seq(dim(correlation_matrix)[2]), function(i){
        removed_na <- na.omit(correlation_matrix[, i])
        # simes_pval <- ifelse(length(removed_na) > 0, simes.test(removed_na), NA)
        if(!(i %in% c(3, 5, 7))){
            pval <- ifelse(length(removed_na) > 0, simes.partial(2, removed_na), NA)
        } else {
            pval <- ifelse(length(removed_na) > 0, extract_binom_pval(removed_na), NA) 
        }
        return(pval)
    })
    return(aggregated_rslt)
}

paralleled.perm.cf <- function(covariates, 
                               Y, treatment, 
                               Y.hat, W.hat, 
                               is_tuned, 
                               num_trees, 
                               mun_perm, 
                               seed = NULL,       # it's seed 
                               is_save = T,
                               file_prefix = NULL, 
                               cluster_id = NULL){

    col_names <- c('obs.Y', 'obs.Y.hat', 'obs.treatment', 'obs.W.hat', 'perm.Y.hat', 'perm.W.hat') 
    cf.estimator <- ifelse(is_tuned, cf.tuned, cf)
    
    perm.risk.mat <- foreach(i = seq(mun_perm), .combine = 'rbind', .options.multicore = list(set.seed = T)) %dopar%  {

        # if cluster_id is provided, we do a stratified permutation to keep cluter heterogeneity present
        if(is.null(cluster_id)){
            samp <- sample(dim(covariates)[1], dim(covariates)[1], replace = F)
        }else{
            samp <- stratified.permutation(cluster_id)
        }

        sampled.covariates <- covariates[samp, ]

        # changed this line 
        perm.tau.forest <- cf.estimator(sampled.covariates, 
                                        Y, treatment, 
                                        Y.hat = Y.hat, 
                                        W.hat = W.hat, 
                                        seed = seed,
                                        num_trees = num_trees)
        if(is_save){
            rslt <- cbind(Y, Y.hat, treatment, W.hat, perm.tau.forest[["Y.hat"]], perm.tau.forest[["W.hat"]])
            ord <- sort(samp, index.return = T)
            write.csv(rslt[ord$ix, ], file = paste0(file_prefix, project, '_repeat_', i,'.csv'), row.names = F, col.names = col_names, quote = F)
        }
        perm.var <- var(perm.tau.forest$predictions)
        perm.risk <- assess.explained.tau.fixed.YW.risk(perm.tau.forest, Y, Y.hat, treatment, W.hat)
        c(perm.var, perm.risk)
    }
    return(perm.risk.mat)
}

permutate.covariates.testing <- function(covariates, Y, 
                                         Y.hat, 
                                         treatment, 
                                         W.hat, 
                                         fixed.YW.tau.risk, 
                                         tau.var, 
                                         is_tuned = F, 
                                         is_save = T, 
                                         file_prefix = NULL, 
                                         cluster_id = NULL, 
                                         num_trees = 1000, 
                                         num.strap = 500, 
                                         seed = NULL){

    # permute covariates with Y, W, Y.hat, and W.hat, to investigate whether covariates contribute to HTE
    # Note: if is_save is T, then file_prefix must be provided. 
    # file_prefix is a prefix for saving, since several results will be saved. 
    # return tau.risk for each permutation
    
    if(is_save && is.null(file_prefix)) stop('if is_save is TRUE, then file_prefix must be provided.')
    
    perm.risk.mat <- paralleled.perm.cf(covariates, 
                                        Y, treatment, 
                                        Y.hat, W.hat, 
                                        is_tuned, 
                                        num_trees, 
                                        num.strap, 
                                        seed = seed, 
                                        is_save = T, 
                                        file_prefix = file_prefix,
                                        cluster_id = NULL)
    if(is_save){    
        observed_risk <- c(tau.var, fixed.YW.tau.risk)
        write.csv(perm.risk.mat, file = paste0(file_prefix, project,  '_fixed_YW_permutation_risk_result.csv'), quote = F, row_names = F)
        write_csv(observed_risk, file = paste0(file_prefix, project, '_fixed_YW_observed_risk_result.csv'), quote = F, row.names = F)
    }
    
    permute.var.pval = mean(perm.risk.mat[,1] >  tau.var)
    fixed.YW.permutation.pval = mean(perm.risk.mat[,2] > fixed.YW.tau.risk)
    
    return(c(permute.var.pval, fixed.YW.permutation.pval))
}

adaptive.permutate.covariates.testing <- function(X, Y, Y.hat, 
                                                 W, W.hat, 
                                                 fixed.YW.tau.risk, 
                                                 tau.var, 
                                                 is_tuned = FALSE, 
                                                 cluster_id = NULL, 
                                                 num_trees = 1000, 
                                                 num.strap = 500, 
                                                 min_perm = 10 ,
                                                 perm_block = 5,
                                                 conf.level = 0.99,
                                                 BinomCI_method = "wilson",
                                                 is_save = T, 
                                                 file_prefix = NULL,
                                                 seed = NULL){
                                                                                                 
    # permute X to investigate whether X contribute to HTE
    # return tau.risk for each permutation
    old.perm.risk.mat <- NULL; stopping <- FALSE; perm_num <- min_perm; idx <- 0

    while(!stopping){

        if(perm_num >= num.strap) break

        perm.risk.mat <- paralleled.perm.cf(X, Y, W, 
                                            Y.hat, W.hat, 
                                            is_tuned, 
                                            num_trees, 
                                            min_perm, 
                                            seed = NULL,   # (seed + perm_num),  
                                            is_save = T, 
                                            file_prefix = paste0(file_prefix, '.', idx), 
                                            cluster_id = NULL)

        if(is.null(old.perm.risk.mat)){
            old.perm.risk.mat <- perm.risk.mat
        } else {
            old.perm.risk.mat <- rbind(old.perm.risk.mat, perm.risk.mat)
        }

        perm_num <- nrow(old.perm.risk.mat)
        
        lower.var.pval <- get.lower.BinomCI(old.perm.risk.mat[,1], tau.var, perm_num, conf.level, BinomCI_method)
        lower.fixed.YW.pval <- get.lower.BinomCI(old.perm.risk.mat[,2], fixed.YW.tau.risk, perm_num, conf.level, BinomCI_method)
        
        idx <- idx + 1
        stopping <- ifelse(early.stop(c(lower.var.pval, lower.fixed.YW.pval)), TRUE, FALSE)
    } # end of while
    
    permute.var.pval <- mean(old.perm.risk.mat[,1] >  tau.var)
    fixed.YW.permute.pval <- mean(old.perm.risk.mat[,2] > fixed.YW.tau.risk)
    
    if(is_save){
        # observed_risk <- c(tau.var, tau.risk, fixed.YW.tau.risk)
        write.csv(old.perm.risk.mat, file = paste0(file_prefix, project, '_permutation_risk_fixed_YW_result.csv'), quote = F, row.names = F)
        # print('b')
        # write.csv(observed_risk, file = paste0(file_prefix, '.observed.risk.result.csv'), quote = F, row.names = F)
    }
    return(c(permute.var.pval, fixed.YW.permute.pval))
}

fisher.exact.test <- function(x, y, sig = 0.1){
    crosstabular  <- table(as.factor(x > sig), as.factor(y > sig))
    unname(crosstabular)
    if(dim(crosstabular)[1] == 2 && dim(crosstabular)[2] == 2){
        fitted <- fisher.test(crosstabular, alternative = 'greater')
        return(fitted$p.value)
    }else{
        return(NA)
    }
}

median.t.test <- function(x, y){
    x.median <- median(x)
    y.median <- median(y)

    if(sum(x <= y.median) < 10 | sum(x > y.median) < 10){
        x.p.value = NA
    }else{
        x.fit <- t.test(x[x <= y.median], x[x > y.median], alternative = 'less')
        x.p.value <- x.fit$p.value
    }

    if(sum(y <= x.median) < 10 | sum(y > x.median) < 10){
        y.p.value = NA
    }else{
        y.fit <- t.test(y[y <= x.median], y[y > x.median], alternative = 'less')
        y.p.value <- y.fit$p.value
    }
    return(c(x.p.value, y.p.value))
}

quantile.t.test <- function(x, y, quantile = 0.5){
    x.quantile <- quantile(x, quantile)
    y.quantile <- quantile(y, quantile)

    if(sum(x <= y.quantile) < 5 | sum(x > y.quantile) < 5){
        x.p.value = NA
    }else{
        x.fit <- t.test(x[x <= y.quantile], x[x > y.quantile], alternative = 'less', na.rm = T)
        x.p.value <- x.fit$p.value
    }

    if(sum(y <= x.quantile) < 5 | sum(y > x.quantile) < 5){
        y.p.value = NA
    }else{
        y.fit <- t.test(y[y <= x.quantile], y[y > x.quantile], alternative = 'less', na.rm = T)
        y.p.value <- y.fit$p.value
    }
    return(c(x.p.value, y.p.value))
}

permutated.pval <- function(Y1, Y2, no.simulations = 10000){
    # test fitting with permutating one outcome variable, and a pval will be returned
    # set.seed(0)
    risk <- mean((Y1 - Y2)^2)

    permutated.risk <- sapply(seq(no.simulations), function(x){
        mean((Y1[sample(length(Y1))] - Y2)^2)
    })
    return(mean(permutated.risk < risk))
}

correlation_test <- function(x, y, methods, alt = 'greater') {

    fitted_res <- sapply(methods, function(mthd) {
        fitted <- cor.test(x, y, method = mthd, alternative = alt, use = "na.or.complete")
        c(fitted$estimate, fitted$p.value)
    })
    fitted_res
}

compute_stats <- function(tau.pred){
    # compute zval, pval and ajusted.p
    tau.zval <- tau.pred$predictions/sqrt(tau.pred$variance.estimates)
    tau.pval <- pnorm(abs(tau.zval), lower.tail = FALSE) * 2
    tau.p.adjust <- p.adjust(tau.pval, method = 'BH')
    stats <- cbind(tau.pred$predictions, tau.zval, tau.pval, tau.p.adjust)
    return(stats)
}
