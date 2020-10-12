format_tcga_patient <- function(pat_ls) {
    tmp = strsplit(pat_ls, "-")
    tmp = unlist(lapply(tmp, function(x) paste(x[[1]], x[[2]], x[[3]], sep = "-")))
    return(unlist(tmp))
}


format_combat <- function(exp_matrix, sample = "tumor") {
    if (sample == "tumor") {
        TP = grepl("-01[[:upper:]]-", exp_matrix$V1)
        exp_matrix = exp_matrix[TP,]
    } else {
        NT = grepl("-01[[:upper:]]-", exp_matrix$V1)
        exp_matrix = exp_matrix[NT,]
    }
    donorId = format_tcga_patient(exp_matrix$V1)
    exp_matrix = cbind(donorId, exp_matrix)
    return(exp_matrix)
}