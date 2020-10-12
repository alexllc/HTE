#' Tutorial
#' from a given casual forest object

#' Step 1: Load the function
source("Code/function.r")

#' Step 2: Provide the forest object file and the training dataset
#' If possible, please make sure that the data and the code are stored
#' in the same folder, otherwise, you may have to modify the regular expression
#' in the get_reprtree function.
forestfile <- "KIRC_ENSG00000066279_forest.RData"
trainingset <- "KIRC_wds.csv"

#' Step 3: Call the get_reprtree function
#' Noticed that this function will return the index of the representative tree
#' to you. You can obtain the tree by command
#' "get_tree(tau.forest, index_tree)"
#' Four distance types were implemented, including d0, d1, d1star, d2
#' Recommend d1star (default).
repretree_index <- get_reprtree(forestfile, trainingset,
                                10, distance_type = "d1star") # Use 10 cores.

#' Do whatever you want with the representative tree index.
