## Main script for survival and ATE analysis of BN genes


library(readxl)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(Homo.sapiens)
library(GeneOverlap)
library(dplyr)
library(NNMIS)
library(grf)
library(data.table)
library(readxl)
library(dplyr)
library(tidyr)

data(GeneOverlap)

source("../HTE/NNMIS_survival_imputation.R")
source("../HTE/HTE_validation_functions.R")

cancer_list = c("PRAD")
# cancer_list = c("BRCA", "PRAD")

for (cancer_abbrv in cancer_list) {
    BN_file_ls = list.files(paste0("./test_files/", cancer_abbrv, "/"))
    for (BN_file in BN_file_ls) {
        message(paste0(c('##', rep('=', 20), ' Processing file: ', BN_file, rep('=', 20)), collapse = ''))
        # Import Bayesian Netowrk selected results
        test_genes = read.csv(paste0("./test_files/", cancer_abbrv, "/", BN_file))
        test_ensembl = test_genes$V1
        test_ensembl = strsplit(test_ensembl, "\\.")
        test_ensembl = sapply(test_ensembl, function(x) {x[[1]][1]})

        # import background gene table for contingency table (file names not hard coded)
        bg_list = list.files(paste0("./bg_gene/", cancer_abbrv))
        whole_blood_index = grep("Whole_Blood", bg_list)
        if(grepl("Whole_Blood", BN_file)) {
            bg_gene = read.csv(paste0("./bg_gene/", cancer_abbrv, "/", bg_list[whole_blood_index]))
        } else {
            bg_gene = read.csv(paste0("./bg_gene/", cancer_abbrv, "/", bg_list[-whole_blood_index]))
        }

        ##### (1) DE GENES SIMILARITIES
        message("##### (1) Comparing with differentially expressed genes.")
        ## 1. Import DE genes result
        DEres = read.csv(paste0("/home/alex/project/HTE/wd/expression_HTE/tables/", cancer_abbrv, "_DEGtable.csv"))
        DE_genes = DEres$X

        # Specify total number of backgroung genes for the contingency tables later
        bg_size = dim(bg_gene)[1]+dim(DEres)[1]

        ## 2. Overlap test
        go.obj <- newGeneOverlap(test_ensembl, DE_genes, genome.size=bg_size) # Background size being total # fo genes tested by 
        go.obj <- testGeneOverlap(go.obj)
        print(go.obj)


        ##### (2) DRIVER GENES SIMILARITIES from Davoli et al. 2013 Table S7. Manually Curated List of Genes Predicted by TUSON Explorer
        message("##### (2) Comparing with pan-cancer driver genes.")
        ## 1. Import driver genes from xlsx
        TSG <- read_excel("davoli_driver_genes.xlsx", sheet = "Table S7A TSGs", skip = 2)
        TSG <- head(TSG[order(TSG$TUSON_q_value_TSG, decreasing= F),], n = 300)
        TSG <- TSG$Gene
        OG <- read_excel("davoli_driver_genes.xlsx", sheet = "Table S7B OGs", skip = 2)
        OG <- head(OG[order(OG$TUSON_q_value_OG, decreasing= F),], n = 250)
        OG <- OG$Gene
        drivers <- c(TSG, OG)
        txnames = AnnotationDbi::select(Homo.sapiens, keys = unique(drivers), columns = "ENSEMBL", keytype = "SYMBOL", multiVals = "CharacterList")
        driver_ensembl = txnames$ENSEMBL

        ## 2. Gene list overlap test
        go.obj <- newGeneOverlap(test_ensembl, driver_ensembl, genome.size=bg_size)
        go.obj <- testGeneOverlap(go.obj)
        print(go.obj)


        ##### Survival advantage for each gene
        message("##### (3) Performing survival analysis.")
        ## 1. Survival time imputation wiht NNMIS

        # Import pre-made OS/expression dataframe
        wds = fread(paste0("/home/alex/project/HTE/wd/expression_HTE/wds_backup/", cancer_abbrv, "_wds.csv"))
        wds = as.data.frame(wds)
        sub_test_ensembl = test_ensembl[test_ensembl %in% colnames(wds)] # overlapping genes between Bayesian Net genes and TCGA data
        sel_col = c("donorId", sub_test_ensembl)
        wds = dplyr::select(wds, all_of(sel_col))

        # Re-imput outcome using NNMIS
        cdr = read_excel("../TCGA-CDR-SupplementalTableS1.xlsx")
        cdr = dplyr::filter(cdr, type == cancer_abbrv)

        # There must be a complete record of OS time and age for NNMIS
        cdr = dplyr::select(cdr, all_of(c("bcr_patient_barcode", "age_at_initial_pathologic_diagnosis", "gender", "OS.time", "OS", "ajcc_pathologic_tumor_stage", "tumor_status"))) %>% filter(!is.na(OS.time)) %>% filter(!is.na(age_at_initial_pathologic_diagnosis)) %>% as.data.frame()
        colnames(cdr)[colnames(cdr) == "bcr_patient_barcode"] = "donorId"

        # The clinical input dataframe for building the survival object must have the same dimension as the dataframe (combined) later used to fit into the Cox PH model.
        common_patients = intersect(cdr$donorId, wds$donorId)
        cdr = cdr[cdr$donorId %in% common_patients,] 

        # Build survival object and impute missing values for auxillary varaiables and survival time
        impute_obj = impute_with_NNMIS(cdr, only_export_obj = TRUE)
        cdr = impute_with_NNMIS(cdr)

        ## 2. Survival analysis
        combined_os_exp = inner_join(cdr, wds, by = "donorId")
    
        # Fit whole dataset into a Cox PH model to check if the genes in question (from BN) significantly affect survival time as a covariate
        attach(combined_os_exp)
        coxest = coxph.pool(obj = impute_obj, time = outcome, status = OS, Z = dplyr::select(combined_os_exp, all_of(sub_test_ensembl)), forceNumeric = FALSE, setRef = NULL)
        message(paste0("Proportion of genes with significant survival prediction ability: ", table(coxest$Pvalue > 0.05)[2]/sum(table(coxest$Pvalue > 0.05)[1], table(coxest$Pvalue > 0.05)[2])))
        detach(combined_os_exp)

        # Combine survival results with test gene ensembl ID
        coxest$var = rownames(coxest) # X represents the effect of stage
        OS_res = test_genes
        OS_res$V1 = test_ensembl
        OS_res = left_join(OS_res, coxest, by = c("V1" = "var"))


        # Only keep cases with complete expression profiles
        # Removing tumor status f we already have stage information as it compromises number of complete cases
        if("ajcc_pathologic_tumor_stage" %in% colnames(combined_os_exp)) {
            combined_os_exp$tumor_status = NULL
        }
        combined_os_exp = combined_os_exp[complete.cases(combined_os_exp),]

        ##### HTE average tx effects
        message("##### (4) Average treatment effect analysis using GRF.")
        ATE_list = data.frame(estimate = NULL, SD = NULL)

        for (gene in sub_test_ensembl) {

            message(paste0("## Processing: ", gene))
            W = as.numeric(combined_os_exp[,gene] > quantile(combined_os_exp[,gene], 0.75)) # using blanket UQ as tx group as most genes from Bayesian net may not be in DE gene list  
            X = dplyr::select(combined_os_exp, -c("donorId", "OS", "OS.time", "outcome"))
            X$gender = as.numeric(as.factor(X$gender))
            Y = combined_os_exp$outcome

            # Build CF (no pre-fitting)
            CF = causal_forest(X,Y,W,
                                Y.hat = NULL,
                                W.hat = NULL,
                                num.trees = 5000,
                                sample.weights = NULL,
                                clusters = NULL,
                                equalize.cluster.weights = FALSE,
                                sample.fraction = 0.5, #When confidence intervals are enabled, the sampling fraction must be less than 0.5.
                                mtry = min(ceiling(sqrt(ncol(X)) + 20), 
                                ncol(X)),
                                min.node.size = 5,
                                honesty = TRUE,
                                honesty.fraction = 0.5,
                                honesty.prune.leaves = TRUE,
                                alpha = 0.05,
                                imbalance.penalty = 0,
                                stabilize.splits = TRUE,
                                ci.group.size = 2,
                                tune.parameters = "none",
                                tune.num.trees = 200,
                                tune.num.reps = 50,
                                tune.num.draws = 1000,
                                compute.oob.predictions = TRUE,
                                orthog.boosting = FALSE,
                                num.threads = 60,
                                seed = 2020
            )
            tau_pred = predict(CF, estimate.variance = TRUE)
            tau_stats = compute_stats(tau_pred)

            # Summarize the p values of all cases with partial simes test (at least 10% of the cases need to be significant in order for this value to be significant)
            simes_pval_ten_percent = simes.partial(floor(dim(tau_stats)[1] * 0.1), tau_stats$tau.pval)
            ATE = c(average_treatment_effect(CF), simes_pval_ten_percent)
            ATE_list = rbind(ATE_list, ATE)
        }

        ATE_res = as.data.frame(cbind(sub_test_ensembl, ATE_list))
        colnames(ATE_res) = c("gene", "ATE_estimate", "ATE_SE", "simes_pval_10_percent")

        all_res = left_join(OS_res, ATE_res, by = c("V1" = "gene"))

        # Convert ensembl versioned genes to gene names
        txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
        txnames =try(AnnotationDbi::select(Homo.sapiens, keys = unique(test_ensembl), columns = "SYMBOL", keytype = "ENSEMBL", multiVals = "CharacterList"))
        all_res = left_join(all_res, txnames, by = c("V1"="ENSEMBL"))

        write.csv(all_res, paste0("./results/", BN_file, "_survival&ATE_res.csv"), row.names = F)
    }
}