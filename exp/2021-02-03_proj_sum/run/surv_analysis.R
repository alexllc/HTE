# Plot basic survival curve with clinical factors for presentations
# tutorial from: https://rviews.rstudio.com/2017/09/25/survival-analysis-with-r/

library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)

# load required data from any HTE R scripts
clin_indexed <- inner_join(clin_indexed, hot_drug, by = "bcr_patient_barcode")


# survival KM object
km <- with(clin_indexed, Surv(OS_time, vital_status))
km_fit <- survfit(km ~ 1, data=clin_indexed)
summary(km_fit, times = c(1,30,60,90*(1:10)))

# KM curve by treatment
km_trt_fit <- survfit(km ~ CETUXIMAB, data=clin_indexed)
summary(km_trt_fit)
summary(km_trt_fit)$table

# Code copied from https://rpkgs.datanovia.com/survminer/
pdf("plot_CETUXIMAB.pdf")
ggsurvplot(
  km_trt_fit,
  data = clin_indexed,
  size = 1,                 # change line size
  palette =
    c("#E7B800", "#2E9FDF"),# custom color palettes
  conf.int = TRUE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  # legend.labs =
  #  c("control", "cetuximab"),    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)
dev.off()

plot(km_trt_fit)