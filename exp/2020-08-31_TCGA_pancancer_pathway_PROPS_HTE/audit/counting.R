library(dplyr)

file_list <- list.files(".")
file_list <- file_list[!grepl(".csv", file_list)]

freq_res <- data.frame(cancer = character(),
                       no_of_tx = numeric(),
                       no_of_perm = numeric(),
                       no_of_corr = numeric(),
                       tc_mean = numeric(),
                       tc_diff = numeric())

for (f in file_list) {
    # perm_res <- read.csv(paste0("./", f, "/PROPS_", f, "_expression_permutate_testing_result.csv"))
    # corr_res <- read.csv(paste0("./", f, "/PROPS_", f, "_expression_correlation_test_result.csv"))
    # tc_res <- read.csv(paste0("./", f, "/PROPS_", f, "_expression_calibration_result.csv"))

    perm_res <- read.csv(paste0("./", f, "/", f, "_expression_permutate_testing_result.csv"))
    corr_res <- read.csv(paste0("./", f, "/", f, "_expression_correlation_test_result.csv"))
    tc_res <- read.csv(paste0("./", f, "/", f, "_expression_calibration_result.csv"))

    # print("Tx analyzied: ")
    no_tx <- dim(corr_res)[1]

    # print("Perm both sig:")
    no_perm <- dim(filter(perm_res, var.pval < 0.05 & fixed.YW.risk.pval < 0.05))[1]

    # print("Corr all sig:")
    no_corr <- dim(filter(corr_res, pearson.pvalue < 0.05 & kendall.pvalue < 0.05 & spearman.pvalue < 0.05))[1]

    no_tc_mean <- dim(filter(tc_res, mean.pred.pval < 0.05))[1]
    no_tc_diff <- dim(filter(tc_res, diff.pred.pval < 0.05))[1]

    row_ls <- do.call("c", c(list(f), as.list(c(no_tx, no_perm, no_corr, no_tc_mean, no_tc_diff))))

    freq_res <- rbind(freq_res, row_ls)
}

colnames(freq_res) <- c("cancer", "no_tx", "no_perm", "no_corr", "no_tc_mean", "no_tc_diff")
write.csv(freq_res, "count_tx.csv", row.names = FALSE)