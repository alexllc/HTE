library(VSURF)

load("covariates_n_gp_for_varselRF.Rdata")
varsel_ranger = VSURF(X, Y, RFimplem = "ranger", parallel = TRUE, verbose = TRUE)
save(varsel_ranger, file = "HNSC_exp_VSURF.RData")