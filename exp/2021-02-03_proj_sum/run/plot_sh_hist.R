# script to visualize the distributions of tau values depending on the overall and correlation significance

# 1. tau p-values significant, correlation significant

# 2. tau p-values not significant, correlation significant

# 3. tau p-values significant, correlation not significant

# 4. tau p-values not significant, correlation not significant

library(dplyr)
library(ggplot2)
library(cowplot)

res_path <- "/home/alex/project/HTE/exp/2020-06-20_TCGA_pancancer_DEA_indicated_tx_dirct_.25_gp_as_tx_exp_HTE/res/BRCA/"

cor_res <- read.csv(paste0(res_path, "BRCA_expression_correlation_test_result.csv"))

case1 <- filter(cor_res, simes.pval < 0.05 & partial.simes.pval < 0.05 & pearson.pvalue < 0.05 & kendall.pvalue < 0.05 & spearman.pvalue < 0.05)

case2 <- filter(cor_res, simes.pval > 0.05 & partial.simes.pval > 0.05 & pearson.pvalue < 0.05 & kendall.pvalue < 0.05 & spearman.pvalue < 0.05)

case3 <- filter(cor_res, simes.pval < 0.05 & partial.simes.pval < 0.05 & pearson.pvalue > 0.05 & kendall.pvalue > 0.05 & spearman.pvalue > 0.05)

case4 <- filter(cor_res, simes.pval > 0.05 & partial.simes.pval > 0.05 & pearson.pvalue > 0.05 & kendall.pvalue > 0.05 & spearman.pvalue > 0.05)

case_ls <- c(case1, case2, case3, case4)

for (case_num in 1:4) {
    
    case <- get(paste0("case", case_num))
    
    gp_a_mat <- matrix(1, nrow = 1076, ncol = length(case$gene))
    gp_b_mat <- matrix(1, nrow = 1076, ncol = length(case$gene))

    for ( k in 1:length(case$gene) ) {
        gene = case$gene[k]
        rep_a <- matrix(1, nrow = 1076, ncol = 10)
        rep_b <- matrix(1, nrow = 1076, ncol = 10)
        for (i in 1:10) {
            sh_a <- read.csv(paste0(res_path, "BRCA_", gene, "_observation_", i, "_result_a.csv"))
            sh_b <- read.csv(paste0(res_path, "BRCA_", gene, "_observation_", i, "_result_b.csv"))
            rep_a[,i] <- sh_a$V1
            rep_b[,i] <- sh_b$V1
        }
        a_mean <- rowMeans(rep_a)
        b_mean <- rowMeans(rep_b)

        gp_a_mat[,k] <- a_mean
        gp_b_mat[,k] <- b_mean
    }

    case_a_mean <- rowMeans(gp_a_mat)
    case_b_mean <- rowMeans(gp_b_mat)

    p1 <- ggplot(as.data.frame(case_a_mean), aes(x = case_a_mean)) + geom_histogram(binwidth = 0.0001, fill="#69b3a2", color="#e9ecef", alpha=0.9) + ggtitle( "First half tau values")
    p2 <- ggplot(as.data.frame(case_b_mean), aes(x = case_b_mean)) + geom_histogram(binwidth = 0.0001, fill="#69b3a2", color="#e9ecef", alpha=0.9) + ggtitle( "Second half tau values")

    pdf(paste0("case_", case_num, "plot.pdf"), width = 8, height = 4)
    print(plot_grid(p1, p2, align = "h"))
    dev.off()

}