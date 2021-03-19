#********************************* TMP *****************************************
exp = "/home/alex/project/HTE/exp/2020-06-20_TCGA_pancancer_DEA_indicated_tx_dirct_.25_gp_as_tx_exp_HTE"
cancer = "BRCA"
perm <- read.csv(paste0(exp, "/res/", "BRCA", "/", "BRCA", "_","expression" , "_permutate_testing_result.csv"))
corr <- read.csv(paste0(exp, "/res/", cancer, "/", cancer, "_", "expression", "_correlation_test_result.csv"))

rep_agg <- data.frame()

for (j in 1:10) {
  sh_a <- read.csv(paste0(exp, "/res/", "/", cancer, "/", cancer, "_", gene_name, "_observation_", j, "_result_a.csv"))
  sh_b <- read.csv(paste0(exp, "/res/", "/", cancer, "/", cancer, "_", gene_name, "_observation_", j, "_result_b.csv"))

  rep_agg <- rbind(rep_agg, sh_a)
  rep_agg <- rbind(rep_agg, sh_b)
}

#****

      rep_tau_psimes <- NULL
      for (j in 1:10) {
        sh_a <- read.csv(paste0(exp, "/res/", "/", cancer, "/", cancer, "_", gene_name, "_observation_", j, "_result_a.csv"))
        sh_b <- read.csv(paste0(exp, "/res/", "/", cancer, "/", cancer, "_", gene_name, "_observation_", j, "_result_b.csv"))

        simes_a <- simes.test(sh_a$V3)
        simes_b <- simes.test(sh_b$V3)

        # FIXME Use Stouffer/Fischer or Fisher's combined probability test for patient-wise independence
        rep_simes <- c(rep_simes, simes_a)
        rep_simes <- c(rep_simes, simes_b)
      }
        
      # For each gene's simes p values, recalculate partial simes for each threshold
      rep_psimes <- get_partial_simes(rep_simes, n = length(rep_simes), u = psimes_u_ls)
      rep_tau_psimes <- rbind(rep_tau_psimes, rep_psimes)
      
      # for each repeat's tau stats, you want to recalculate 
      for (pcols in )

    }
    # write summarized tables to excel sheets
    write.xlsx(perm, file = "summary_output.xlsx", sheetName = paste0("DE", opt$direction, "_permutation"), append = FALSE)
    write.xlsx(perm, file = "summary_output.xlsx", sheetName = paste0("DE", opt$direction, "_correlation"), append = TRUE)
    write.xlsx(tau_tbl, file = "summary_output.xlsx", sheetName = paste0("DE", opt$direction, "_tau_values"), append = TRUE)
