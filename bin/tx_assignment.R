
#' Function for formatting treatment assignment
#' Will require "./dat/tables" to be present
#' 
#' @param binary [logical] whether treatment group assignment should be made binary
#' @param upperQ [logical] if binary = TRUE, whether set upper quantile as the treated W = 1 or untreated W = 0 group.
#' @param thres [numeric] where the cut-off for upperQ should be, default is 0.75.
#' @param treatment [numeric vector] of original per patient entries, could be pharmacological, non-pharm, proxy (gene expression/mutations) etc.
#' @return [numeric vector] formatted treatment vector

assign_tx <- function(binary = TRUE, upperQ = NULL, thres = 0.75, treatment = NULL) {
    if(binary){
        if (!is.null(upperQ) ) {
            if(upperQ) {
                treatment = as.numeric(treatment > quantile(treatment, thres))
            } else {
                treatment  = as.numeric(treatment < quantile(treatment, 1- thres))
            }
        } else {
            treatment = as.numeric(treatment != 0)
        } # only for mutation
    } else {
        treatment = as.numeric(treatment)
    }
    return(treatment)
}