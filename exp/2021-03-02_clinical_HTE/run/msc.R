# Pancancer info

# Check date of expression vs date of drug administration
# > colnames(drug)[grep("date", colnames(drug))]
# [1] "form_completion_date"

library(readxl)
library(RColorBrewer)
library(dplyr)

cancerList <- c("BLCA", "BRCA", "COAD", "GBM", "HNSC", "LUAD", "LUSC", "KIRC", "PRAD", "STAD" ,"UCEC")

fit <- 0

for (cancer_type in cancerList) {
    res <- read_excel(paste0("/home/alex/project/HTE/exp/2021-03-02_clinical_HTE/res_raw_OS_time/", cancer_type, "/", cancer_type, "_grf_res.xlsx"), sheet = "test_calib")
    # message(paste0("No of drugs examined for ", cancer_type))
    # print(dim(res)[1])

    res[["Pr..t."]] %<>% as.numeric

    message(paste0("No of sig models for ", cancer_type))
    print(dim(filter(res, `Pr..t.` < 0.05))[1])
    fit <- fit + dim(filter(res, `Pr..t.` < 0.05))[1]

    # sig ATE

}

for (i in 2:dim(res)[2]) {
    res[[colnames(res)[i]]] %<>% as.numeric
}

# get ATE extract into cleaveland plot

agg_ate <- data.table()
ate <- 0
for (cancer_type in cancerList) {
    res <- read_excel(paste0("/home/alex/project/HTE/exp/2021-03-02_clinical_HTE/res_raw_OS_time/", cancer_type, "/", cancer_type, "_grf_res.xlsx"), sheet = "avg_tx_eff")
    for (i in 2:dim(res)[2]) {
        res[[colnames(res)[i]]] %<>% as.numeric
    }
    ate_drugs <- which(abs(res$estimate) > res$std.err * 2)
    no_of_ate <- length(ate_drugs)
    types <- rep(cancer_type, no_of_ate)
    save_ate <- res[ate_drugs,]
    save_ate <- bind_cols(types, save_ate)
    colnames(save_ate)[1] <- "cancer"
    agg_ate <- bind_rows(agg_ate, save_ate)

    print(no_of_ate)
    ate <- ate + no_of_ate

    message(paste0("No of sig ATE for ", cancer_type))
}


#### PLOT

pdf(file = "ATE_plot.pdf")
ggdotchart(agg_ate, x = "V1", y = "estimate",
           color = "cancer",                                # Color by groups
           palette = brewer.pal(length(unique(agg_ate$cancer)), "Spectral"), # Custom color palette
           sorting = "descending",                       # Sort value in descending order
           rotate = TRUE,                                # Rotate vertically
           dot.size = 3,                                 # Large dot size
           y.text.col = TRUE, 
           font.label =  list(size = 11),                           # Color y text by groups
           ggtheme = theme_pubr()                        # ggplot2 theme
           )+
  theme_cleveland()                                      # Add dashed grids
  dev.off()

df = as.data.frame(agg_ate[order(agg_ate$estimate, decreasing = TRUE),])
pdf(file = "ATE_plot.pdf")
ggplot(
  df, 
  aes(x = estimate, y = reorder(V1, estimate), xmin = estimate - std.err, xmax = estimate + std.err)
  ) +
  geom_point(aes(color = cancer)) +
  geom_errorbarh(aes(color = cancer), height=.2)+
  labs(x = "Tau estimate (days)", y = "Drug type") +

  theme_light()
dev.off()


## PLOT TAU

tau <- read.csv("/home/alex/project/HTE/exp/2021-03-02_clinical_HTE/res_raw_OS_time/LUAD/perm_shc/LUAD_tau_PEMETREXED.csv")
tau <- read.csv("/home/alex/project/HTE/exp/2021-03-02_clinical_HTE/res_raw_OS_time/BRCA/perm_shc/BRCA_tau_ANASTROZOLE.csv")

tau_sub <- filter(tau, tau.p.adjust < 0.05)
tau_sub <- tau_sub[!duplicated(tau_sub$donorId),]
tau_sub$direction <- ifelse(tau_sub$tau.val > 0, "Positive", "Negative")
tau_sub$percentage <- exp(abs(tau_sub$tau.val)) * 100

pdf("tau_zscore.pdf")
ggbarplot(tau_sub, x = "donorId", y = "percentage",
          fill = "direction",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in ascending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Tau Z-score",
          xlab = FALSE,
          legend.title = "MPG Group"
          )
dev.off()

pdf("tau_zscore_dot.pdf")
ggdotchart(tau_sub, x = "donorId", y = "tau.zval",
           color = "direction",                                # Color by groups
           palette = c("#E7B800", "#FC4E07"), # Custom color palette
           sorting = "ascending",                        # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           ggtheme = theme_pubr(),
           xlab = FALSE,
           ylab = "Tau Z score"                        # ggplot2 theme
           ) +
             theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                      axis.ticks.x=element_blank())
dev.off()

## Varimp enrichements plot

varimp <- read.csv("")
class(varimp$importance) <- "numeric"
varimp_gene <- varimp[grep("ENSG*", varimp$feature),]
varimp_clin <- varimp[!grep("ENSG*", varimp$feature),] # only keep genes
# Convert ensemblIDs to ENTERZ gene ID for the query
cutoff <- quantile(varimp_gene$importance, 0.75)
varimp_gene <- filter(varimp_gene, importance > cutoff)
suppressMessages(ensembl2ID <- AnnotationDbi::select(org.Hs.eg.db, keys=varimp$feature,columns=c("ENTREZID"), keytype="ENSEMBL")) 

kk <- enrichKEGG(gene = ensembl2ID$ENTREZID,
                organism = 'hsa',
                pvalueCutoff = 0.05)
# kk <- filter(as.data.frame(kk), p.adjust < 0.05)
# k_id <- paste(kk$ID, collapse = "|")
# k_names <- paste(kk$Description, collapse = "|")
# top_varimp <- rbind(top_varimp, c(tx_gene_name, k_id, k_names))

# plot_kk <- dplyr::select(kk, c("Description", "GeneRatio"))
# tmp <- strsplit(plot_kk$GeneRatio, "/")
# plot_kk$num <- as.numeric(sapply(tmp, function(x) x[[1]]))
# plot_kk$denom <- as.numeric(sapply(tmp, function(x) x[[2]]))
# plot_kk$percent <- plot_kk$num / plot_kk$denom * 100

library(enrichplot)

pdf("KEGG.pdf")
barplot(kk, showCategory=20) +
    ggtitle("Anastrozole varimp") +
    labs(x = "KEGG pathway names", y = "Gene ratio")
dev.off()


## back elim
elim_file <- "/home/alex/project/HTE/exp/2021-03-02_clinical_HTE/res/LUAD/perm_shc/LUAD_back_elim_blp_PEMETREXED.csv"
elim_file <- "/home/alex/project/HTE/exp/2021-03-02_clinical_HTE/res/BRCA/perm_shc/BRCA_back_elim_blp_ANASTROZOLE.csv"
saved_elim <- read.csv(elim_file)


elim_gene <- saved_elim$X[grep("ENSG*", saved_elim$X)]
suppressMessages(ensembl2names <- AnnotationDbi::select(org.Hs.eg.db, keys=elim_gene,columns=c("SYMBOL"), keytype="ENSEMBL"))
named_elim <- left_join(saved_elim, ensembl2names, by = c("X" = "ENSEMBL"))
write.csv(named_elim, file = elim_file, row.names = FALSE)

named_elim$SYMBOL[1] <- "Intercept"
named_elim$direction <- ifelse(named_elim$Estimate > 0, "Positive", "Negative")

pdf("back_elim.pdf")
ggbarplot(named_elim, x = "SYMBOL", y = "Estimate",
          fill = "direction",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "desc",          # Sort the value in descending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Best linear prediction estimate",
          xlab = "Covariates",
          legend.title = "Direction",
          rotate = TRUE,
          ggtheme = theme_minimal()
          )
dev.off()