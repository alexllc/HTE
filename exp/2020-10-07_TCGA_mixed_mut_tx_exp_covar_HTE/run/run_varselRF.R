library(varSelRF)

load("mut_covariates_n_gp_for_varselRF.Rdata")
set.seed(2020)
varsel = varSelRF(X, as.factor(Y))
save(varsel, file = "HNSC_mutation_varsel_results.Rdata")

# Expression covariates
# check 75335
## saved Oct  7 18:53, started Oct  8 07:25

# Mutation covariates
# Check 4012