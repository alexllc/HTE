2020-10-14 17:30

# Meeting with Prof So and Larry

## Objectives:
1. Mutation HTE can't be used as it is (data too sparse and the trees are not reliable)
2. Using NLP or RF to pre-select variables
3. High censoring problem in cancer types requires alternative outcome measures


## Discussion:
1. Use expression covariate matrix in place of mutation matrix
2. Interesting direction but for the next paper or project
   
   Pre-selection inevitably introduce bias. You could do post-selection inference, which is an area of research not a lot of computer scientists are interested in, but we can look into it.

   NLP is useful for validating models, could pair with LASSO, but Larry's current NLP doesn't have a weighting or score yet.

   The weight can be calculated simply by formula $\gamma_j = \frac{1}{\text{cite-count } g_j + 1}$. $g_j$ means gene $j$.
   Noticed that a pseudo-count is added here.  
   Therefore, we have the lasso in **lagrangian form**,  
   $$\min_{\beta\in\mathbb{R}^P}\{\frac{1}{N} \lVert y-X\beta\lVert^2_2\} + \lambda\lVert \gamma\beta \lVert_1$$

3. Try **Time to remission**, **Progression-free survival** or **Overall survival** for the new mutation/expression combination HTE analysis -> for next discussion



## To-do list:

- [x] Assess the amount of censoring across cancer types 
   See `exp/2020-10-07_TCGA_mixed_mut_tx_exp_covar_HTE/dat/proportion_of_censoring.csv`, script under `exp/2020-10-07_TCGA_mixed_mut_tx_exp_covar_HTE/run/pancancer_censoring.R`
- [x] retreive alternative outcome measures & establish appropriate time to failure function
- [x] Perform combined mutation/expression HTE using the three mentioned outcome for further comparison
- [x] Re-do external validation with mut/exp combination HTE
