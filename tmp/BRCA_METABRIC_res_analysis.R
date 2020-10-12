f <- list.files(getwd(), pattern="\\_split_half.csv$")

correlation_test_ret = NULL
corr_test_names = c("gene", "simes.pval", "partial.simes.pval", "pearson.estimate", "pearson.pvalue", "kendall.estimate", "kendall.pvalue", "spearman.estimate", "spearman.pvalue")

for (file in f){
    correlation_matrix = as.matrix(read.csv(file))
    aggregated_corr_rslt <- sapply(seq(dim(correlation_matrix)[2]), aggr_res, est_col_list = c(3, 5, 7), res_mat = correlation_matrix)
    tx = strsplit(file, "_")[[1]][2]
    current_ret <- do.call('c', list(list(tx), as.list(aggregated_corr_rslt[1: 8])))
    current_ret <- rbindlist(list(current_ret))
    correlation_test_ret <- rbind(correlation_test_ret, current_ret)
}

colnames(correlation_test_ret) = corr_test_names

write.csv(correlation_test_ret, paste0(project, '_validation_correlation_test_result.csv'), quote = F, row.names = F)