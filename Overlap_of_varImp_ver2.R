# compare the concordnace of VarImp in two datasets (eg split-half) 
# project/HTE/wd/mut_HTE/result/BRCA_METABRIC_TCGA_mut_freq_validation

# Input matrix 
# 1st col: gene_name /ID
# 2nd col: varImp of first dataset 
# 3rd col: VarImp of 2nd dataset 
library(foreach)
library(TFisher) 
library(metap)
library(parallel)


degree_of_overlap <- function(varimp_mat, cutoff = c(0,0), estimate = F) {
    genels = varimp_mat[,1]
    rank_A = order(abs(varimp_mat[,2]),decreasing = TRUE)
    rank_B = order(abs(varimp_mat[,3]),decreasing = TRUE)
    A = genels[rank_A][1:cutoff[1]]
    B = genels[rank_B][1:cutoff[2]]
    contingen_tbl = fisher.test(matrix(c(length(intersect(A,B)), length(setdiff(A,B)), length(setdiff(B,A)), length(genels) - length(union(A,B))), nrow=2), alternative = "greater")
    if (estimate) {
        return(c(contingen_tbl$p.value, contingen_tbl$estimate))
    } else {
       return(contingen_tbl$p.value)
    }
}

simes.test <- function(x, returnstat=FALSE){
    r=rank(x,  ties.method = "random")
    T=min(length(x)*x/r)
    if (returnstat) c(T,T) else T
}

test_overlap_VarImp <- function (input_matrix,
								no_perm = 5000,
								no_cluster = 10, #no. of clusters for parallel running 
								top_percentile_list = c(0.01, 0.03, 0.05, 0.1, 0.2, 0.3) ) {

    list_quantiles = quantile(1:dim(input_matrix)[1], top_percentile_list )
    list_TopK_genes = round(list_quantiles)   #c(20, 50, 100, 200, 300) 
        
    fisher_test_result = NULL
    fisher_res_topK = NULL
    for (i in 1:length(list_TopK_genes)) {      
        fisher_res_topK[i] = degree_of_overlap(input_matrix, cutoff = as.numeric(rep(list_TopK_genes[i], 2)))
    }

    ## Degree of overlap between genes with varimp > 0
    fisher_res_above0 = degree_of_overlap(input_matrix, cutoff = c(sum(input_matrix[,2] > 0), sum(input_matrix[,3] > 0))) # originally: getting estimate too
    # names(fisher_res_above0) = c("pval", "est")

    ## Combining p values across quantiles and >0
    pall = c(fisher_res_topK, as.numeric(fisher_res_above0[1]))
    obs_fisher =  sumlog(pall)$p
    obs_min = min(pall)
    obs_simes =  simes.test(pall)	
    TAU1 = c(0.05, 0.1, 0.5, 1)
    no_p = length(pall)
    q_omni = stat.soft.omni(p=pall, TAU1=TAU1)
    obs_softomni = p.soft.omni(q=q_omni$omni, n=no_p, TAU1=TAU1, M = NULL)

    ## Permutation tests
    perm_fisher=numeric(no_perm)
    perm_min=numeric(no_perm)
    perm_simes=numeric(no_perm)
    perm_softomni=numeric(no_perm)

    cl <- parallel::makeCluster(no_cluster, type = "FORK")
    doParallel::registerDoParallel(cl)

    all_perm_results = foreach(j = 1:no_perm, .combine = 'rbind', .export=c('degree_of_overlap', 'simes.test'), .packages= c("TFisher", "metap")) %dopar% {
        perm_matrix = cbind(input_matrix[,c(1,2)], sample(input_matrix[,3]))
        perm_fisher_res = NULL

        for (i in 1:length(list_TopK_genes)) {
            perm_fisher_res[i] = degree_of_overlap(perm_matrix, cutoff = as.numeric(rep(list_TopK_genes[i], 2)))
        }

        perm_fisher_res_above0 = degree_of_overlap(perm_matrix, cutoff = c(sum(perm_matrix[,2] > 0), sum(perm_matrix[,3] > 0)))

        pall = c(perm_fisher_res, perm_fisher_res_above0)
        perm_fisher =  sumlog(pall)$p
        perm_min = min(pall)
        perm_simes =  simes.test(pall)	

        TAU1 = c(0.05, 0.1, 0.5, 1)
        no_p = length(pall) 
        q_omni = stat.soft.omni(p=pall, TAU1=TAU1)
        perm_softomni = p.soft.omni(q=q_omni$omni, n=no_p, TAU1=TAU1, M = NULL)
        perm_list = c(perm_fisher, perm_min, perm_simes, perm_softomni)
    }

    parallel::stopCluster(cl)

    perm_p1 = (sum(all_perm_results[,1] < obs_fisher)+1)/(no_perm+1)
    perm_p2 = (sum(all_perm_results[,2]  < obs_min)+1)/(no_perm+1)
    perm_p3 = (sum(all_perm_results[,3]  < obs_simes )+1)/(no_perm+1)
    perm_p4 = (sum(all_perm_results[,4]  < obs_softomni)+1)/(no_perm+1)

    perm_p_res = c(perm_p1, perm_p2, perm_p3, perm_p4) 

    return( list(fisher_res_topK  = fisher_res_topK ,
                top_percentile_list = top_percentile_list ,
                fisher_res_above0 = fisher_res_above0,
                perm_p_combineALL = perm_p_res) ) 
			 
}