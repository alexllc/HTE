# compare the concordnace of VarImp in two datasets (eg split-half) 

# Input matrix 
# 1st col: gene_name /ID
# 2nd col: varImp of first dataset 
# 3rd col: VarImp of 2nd dataset 
library(foreach)
library(TFisher) 
library(metap) 


#example input 
# file1= read.csv("/exeh_3/yuping/demap_expression/BRCA/demap_BRCA_expression/FAM110D_interaction_test_expr.csv")
# file2= read.csv("/exeh_3/yuping/demap_expression/BRCA/demap_BRCA_expression/JPH3_interaction_test_expr.csv")
# input_matrix = merge(file1, file2, by="Ensembl_Gene_ID")
# input_matrix= input_matrix[,c(1,3,8)]


test_overlap_VarImp <- function (input_matrix,
								no_perm = 5000,
								no_cluster = 10, #no. of clusters for parallel running 
								top_percentile_list = c(0.01, 0.03, 0.05, 0.1, 0.2, 0.3) ) {

    #library(energy)
    #library(HHG) 




    ##End of Input_____________________________



    simes.test <- function(x, returnstat=FALSE){
    r=rank(x,  ties.method = "random")
    T=min(length(x)*x/r)
    if (returnstat) c(T,T) else T
    }



    gene_list = input_matrix[,1]
    imp1 = as.numeric( input_matrix[,2] ) 
    imp2 = as.numeric( input_matrix[,3] ) 

    genes_num = length(gene_list) 
    list_quantiles = quantile(1:genes_num, top_percentile_list )
    list_TopK_genes = round(list_quantiles)   #c(20, 50, 100, 200, 300) 

    #add some rank inforamtion

    order_1st = order(abs(imp1),decreasing = TRUE)
    order_2nd = order(abs(imp2),decreasing = TRUE)
        
        #show dataframe in order with a particular column
        # dt_IDA2 = dt_IDA[order(dt_IDA$Rank_ABS_true_IDA),]
        
    fisher_test_result = NULL
    fisher_res_topK = NULL
    for (i in 1:length(list_TopK_genes)) {      
        TopK_genes = list_TopK_genes[i]
        top_gene_1st_list = gene_list[order_1st][1:TopK_genes] # obtain topKgenes in trueIDA list
        top_gene_2nd_list = gene_list[order_2nd][1:TopK_genes] #obtain topK genes in estimateIDA list
        
        #top_gene_1st_list = order_1st[1:TopK_genes] # obtain topKgenes in trueIDA list
        #top_gene_2nd_list = order_2nd[1:TopK_genes] #obtain topK genes in estimateIDA list
        
        intersection = intersect(top_gene_1st_list,top_gene_2nd_list ) # obtain the intersection of these two list, and the series stand the gene serial number
        
        contingency_table_11 = length(intersection)
        contingency_table_12 = TopK_genes-contingency_table_11
        contingency_table_21 = TopK_genes-contingency_table_11
        contingency_table_22 = genes_num-contingency_table_21-contingency_table_11-contingency_table_12
        
        fisher_res_topK[i] = fisher.test(rbind(c(contingency_table_11,contingency_table_12),c(contingency_table_21,contingency_table_22)), alternative="greater")$p.value
        
    }

    #fisher_res_topK 


    top_gene_1st_list = gene_list[imp1>0] # obtain topKgenes in trueIDA list
    top_gene_2nd_list = gene_list[imp2>0] #obtain topK genes in estimateIDA list
        
    intersection = intersect(top_gene_1st_list,top_gene_2nd_list ) # obtain the intersection of these two list, and the series stand the gene serial number
        
    contingency_table_11 = length(intersection)
    contingency_table_12 =  sum(imp2>0) - contingency_table_11
    contingency_table_21 =  sum(imp1>0) - contingency_table_11
    contingency_table_22 = genes_num -contingency_table_21-contingency_table_11-contingency_table_12
        
    fisher_res_aboveZero.obj = fisher.test(rbind(c(contingency_table_11,contingency_table_12),c(contingency_table_21,contingency_table_22)), alternative="greater")
    fisher_res_aboveZero <-  c(fisher_res_aboveZero.obj$p.value, fisher_res_aboveZero.obj$estimate) 

    #__________________________________________________________
    #  combine across p-values at differnet quantiles 
    #____________________________________________________________

    pall = c(fisher_res_topK, fisher_res_aboveZero) 
    obs_fisher =  sumlog(pall)$p
    obs_min = min(pall)
    obs_simes =  simes.test(pall)	

    TAU1 = c(0.05, 0.1, 0.5, 1)
    no_p = length(pall) 
    q_omni = stat.soft.omni(p=pall, TAU1=TAU1)
    obs_softomni = p.soft.omni(q=q_omni$omni, n=no_p, TAU1=TAU1, M = NULL)







    perm_fisher=numeric(no_perm)
    perm_min=numeric(no_perm)
    perm_simes=numeric(no_perm)
    perm_softomni=numeric(no_perm)

    cl <- parallel::makeCluster(no_cluster)
    doParallel::registerDoParallel(cl)

    t1 = proc.time()

    all_perm_results = foreach(j = 1:no_perm, .combine = 'rbind') %dopar% {
        library(TFisher) 
        library(metap) 

        #_______________________________________________
        # permutation to determine signicance
        #_______________________________________________


        #to do: permutation: imp1 <- sample(imp1) 
        # combine pval across thresholds by fisher, min, omni



        gene_list = input_matrix[,1]

        #_______________________________________________
        # note the permutation here
        #______________________________________________
        imp1 = sample(  as.numeric( input_matrix[,2] )   ) 
        imp2 = as.numeric( input_matrix[,3] ) 

        genes_num = length(gene_list) 
        list_quantiles = quantile(1:genes_num, top_percentile_list )
        list_TopK_genes = round(list_quantiles)   #c(20, 50, 100, 200, 300) 

        #add some rank inforamtion

        order_1st = order(abs(imp1),decreasing = TRUE)
        order_2nd = order(abs(imp2),decreasing = TRUE)
            
            #show dataframe in order with a particular column
            # dt_IDA2 = dt_IDA[order(dt_IDA$Rank_ABS_true_IDA),]
            
        fisher_test_result = NULL

        for (i in 1:length(list_TopK_genes)) {      
            TopK_genes = list_TopK_genes[i]
            top_gene_1st_list = gene_list[order_1st][1:TopK_genes] # obtain topKgenes in trueIDA list
            top_gene_2nd_list = gene_list[order_2nd][1:TopK_genes] #obtain topK genes in estimateIDA list
            
            #top_gene_1st_list = order_1st[1:TopK_genes] # obtain topKgenes in trueIDA list
            #top_gene_2nd_list = order_2nd[1:TopK_genes] #obtain topK genes in estimateIDA list
            
            intersection = intersect(top_gene_1st_list,top_gene_2nd_list ) # obtain the intersection of these two list, and the series stand the gene serial number
            
            contingency_table_11 = length(intersection)
            contingency_table_12 = TopK_genes-contingency_table_11
            contingency_table_21 = TopK_genes-contingency_table_11
            contingency_table_22 = genes_num-contingency_table_21-contingency_table_11-contingency_table_12
            
            fisher_test_result[i] = fisher.test(rbind(c(contingency_table_11,contingency_table_12),c(contingency_table_21,contingency_table_22)), alternative="greater")$p.value
            
        }

        #fisher_test_result


        top_gene_1st_list = gene_list[imp1>0] # obtain topKgenes in trueIDA list
        top_gene_2nd_list = gene_list[imp2>0] #obtain topK genes in estimateIDA list
            
        intersection = intersect(top_gene_1st_list,top_gene_2nd_list ) # obtain the intersection of these two list, and the series stand the gene serial number
            
        contingency_table_11 = length(intersection)
        contingency_table_12 =  sum(imp2>0) - contingency_table_11
        contingency_table_21 =  sum(imp1>0) - contingency_table_11
        contingency_table_22 = genes_num -contingency_table_21-contingency_table_11-contingency_table_12
            
        fisher_test_cutoffzero = fisher.test(rbind(c(contingency_table_11,contingency_table_12),c(contingency_table_21,contingency_table_22)), alternative="greater")$p.value
        #fisher_test_cutoffzero   

        #__________________________________________________________
        #  combine across p-values at differnet quantiles 
        #____________________________________________________________

        pall = c(fisher_test_result, fisher_test_cutoffzero) 
        perm_fisher =  sumlog(pall)$p
        perm_min = min(pall)
        perm_simes =  simes.test(pall)	

        TAU1 = c(0.05, 0.1, 0.5, 1)
        no_p = length(pall) 
        q_omni = stat.soft.omni(p=pall, TAU1=TAU1)
        perm_softomni = p.soft.omni(q=q_omni$omni, n=no_p, TAU1=TAU1, M = NULL)

        ##results of fisher p, min p etc. under the permutation null distribution
        perm_list = c(perm_fisher, perm_min, perm_simes, perm_softomni) 
        #perm_list  

    }

    #save(all_perm_results, file="/exeh_4/sohc/heterogenous_treatment_effects/Dependency_map_validation/BRCA_permuted_results.Rdata") 
    t2 = proc.time()
    t2-t1
    parallel::stopCluster(cl)



    #______________________________________________
    # final p-values 
    #____________________________________________


    perm_p1 = (sum(all_perm_results[,1] < obs_fisher)+1)/(no_perm+1)
    perm_p2 = (sum(all_perm_results[,2]  < obs_min)+1)/(no_perm+1)
    perm_p3 = (sum(all_perm_results[,3]  < obs_simes )+1)/(no_perm+1)
    perm_p4 = (sum(all_perm_results[,4]  < obs_softomni)+1)/(no_perm+1)

    perm_p_res = c(perm_p1, perm_p2, perm_p3, perm_p4) 


    # perm_p_res 


    return( list(fisher_res_topK  = fisher_res_topK ,
                top_percentile_list = top_percentile_list ,
                fisher_res_aboveZero = fisher_res_aboveZero,
                perm_p_combineALL = perm_p_res) ) 
			 
}





# res = test_overlap_VarImp (input_matrix,
# 								no_perm = 5000,
# 								no_cluster = 10, #no. of clusters for parallel running 
# 								top_percentile_list = c(0.01, 0.03, 0.05, 0.1, 0.2, 0.3) 
								
# 					)
					
# res