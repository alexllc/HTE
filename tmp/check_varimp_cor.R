library(grf)
library(plyr)
library(dplyr)
library(tidyr)

 X.covariates <- as.matrix(dplyr::select(covar_mat, -tx))
 
 forest = causal_forest(X.covariates, Y, treatment, mtry = round(ncol(X.covariates)* 0.01), num.trees = 2000, num.threads=10)

 ate = average_treatment_effect(forest)
 tree.plot = plot(get_tree(forest, 10))
 cat(DiagrammeRsvg::export_svg(tree.plot), file='CD8_gene_plot3.svg')

 covar = variable_importance(forest)
 covar = data.frame(colnames(X.covariates), var)
 covar = var[order(var$var, decreasing=T),]

####################################################################
varimp = read.csv("./result/HNSC/HNSC_varimp_CD8_Tcells.csv")
varimp$variable = as.character(varimp$variable)
tau = read.csv("./expression_HTE/result/BRCA/BRCA_tau_ENSG00000205364.csv")
genevar = varimp[varimp$variable %in% colnames(exp_matrix),]

res = NULL
for (covariate in genevar$variable) {
     covar_value = exp_matrix[,colnames(exp_matrix)==covariate]
     res =cbind(res, cor.test(tau$tau.val, as.numeric(covar_value)[1:1089]))
 }
colnames(res) = genevar$variable
res = t(res)
res = as.data.frame(res)
res$lower.conf = unlist(lapply(res$conf.int, function(x) x[[1]]))
res$upper.conf = unlist(lapply(res$conf.int, function(x) x[[2]]))
res$conf.int = NULL
res = apply(res, 2, function(x) do.call("rbind", x))
res = as.data.frame(apply(res, 2, as.numeric))
res$variable = genevar$variable


varimp_cor = left_join(varimp, res, by = "variable")
write.csv(varimp_cor, "BRCA_MT1M_varimp_cor.csv")