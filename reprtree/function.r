#' Obtain binary strings.
#' @tree casual tree object.
#' @data the data used to train the casual forest.
obtain_bs <- function(tree, data) {
    variable_label <- rep(0, length(colnames(data)))

    nodes  <-  tree$nodes
    is_leaf <- sapply(nodes, function(x) x["is_leaf"])
    ind_leaf <- which(unlist(is_leaf) == FALSE)
    node_split_set <- sapply(nodes, function(x) x["split_variable"])
    node_split_set <- unlist(node_split_set[ind_leaf])
    variable_label[node_split_set] <- 1

    return(variable_label)
}

#' Obtain Hamming distance between two trees. (d0)
#' @binlist1 binary strings corresponding to tree 1.
#' @binlist2 binary strings corresponding to tree 2.
obtain_d0_distance <- function(tree1, tree2, data) {
    binlist1 <- obtain_bs(tree1, data)
    binlist2 <- obtain_bs(tree2, data)

    d_0_t1_t2 <- length(union(setdiff(binlist1, binlist2),
                        setdiff(binlist2, binlist1)))/length(colnames(data))

    return(d_0_t1_t2)
}

#' Obtain d1 distance between two trees.
#' @tree1 tree 1.
#' @tree2 tree 2
obtain_d1_distance <- function(tree1, tree2) {
    # Get the nodes list of tree 1 and tree 2 respectively
    nodes1  <-  tree1$nodes
    nodes2  <-  tree2$nodes

    # Obtain the subjects that are included in both trees
    ovlap_subj <- sort(intersect(tree1$drawn_samples, tree2$drawn_samples))
    total_pair <- as.data.frame(t(combn(ovlap_subj, 2)))

    # Obtain the subject selected by both trees during construction
    tree1_sample <- sort(intersect(tree1$drawn_samples, ovlap_subj))
    tree2_sample <- sort(intersect(tree2$drawn_samples, ovlap_subj))

    # Obtain the pair subject not located in the same node in both tree
    is_leaf <- sapply(nodes1, function(x) x["is_leaf"])
    ind_leaf <- which(unlist(is_leaf) == TRUE)
    cluster_samp <- sapply(nodes1, function(x) x["samples"])
    cluster_samp <- cluster_samp[ind_leaf]
    cluster_samp <- lapply(cluster_samp, function(x) intersect(x, ovlap_subj))

    ind_not_zero <- which(unlist(sapply(cluster_samp,
                                        function(x) length(x))) != 0)
    cluster_samp <- cluster_samp[ind_not_zero]

    cluster_samp_unique <- lapply(cluster_samp, function(x) t(combn(x,2)))
    df_1 <- distinct(as.data.frame(do.call("rbind", cluster_samp_unique)))

    is_leaf <- sapply(nodes2, function(x) x["is_leaf"])
    ind_leaf <- which(unlist(is_leaf) == TRUE)
    cluster_samp <- sapply(nodes2, function(x) x["samples"])
    cluster_samp <- cluster_samp[ind_leaf]
    cluster_samp <- lapply(cluster_samp, function(x) intersect(x, ovlap_subj))

    ind_not_zero <- which(unlist(sapply(cluster_samp,
                                        function(x) length(x))) != 0)
    cluster_samp <- cluster_samp[ind_not_zero]

    cluster_samp_unique <- lapply(cluster_samp, function(x) t(combn(x, 2)))
    df_2 <- distinct(as.data.frame(do.call("rbind", cluster_samp_unique)))

    df_3 <- data.frame(V1 = df_2$V2, V2 = df_2$V1)

    a <- intersect(df_1, df_2)
    b <- intersect(df_1, df_3)
    c <- rbind(a, b)

    numerator <- dim(df_1)[1] + dim(df_2)[1] - 2 * dim(c)[1]
    denominator <- choose(length(ovlap_subj), 2)

    return(numerator/denominator)
}

#' Obtain d1* distance between two trees.
#' @tree1 tree 1.
#' @tree2 tree 2
obtain_d1_star_distance <- function(tree1, tree2, data) {
    # Get the nodes list of tree 1 and tree 2 respectively
    nodes1  <-  tree1$nodes
    nodes2  <-  tree2$nodes

    # Obtain the subjects that are included in both trees
    ovlap_subj <- sort(intersect(tree1$drawn_samples, tree2$drawn_samples))
    total_pair <- as.data.frame(t(combn(ovlap_subj, 2)))

    # Obtain the subject selected by both trees during construction
    tree1_sample <- sort(intersect(tree1$drawn_samples, ovlap_subj))
    tree2_sample <- sort(intersect(tree2$drawn_samples, ovlap_subj))

    # Obtain the pair subject not located in the same node in both tree
    is_leaf <- sapply(nodes1, function(x) x["is_leaf"])
    ind_leaf <- which(unlist(is_leaf) == TRUE)
    cluster_samp <- sapply(nodes1, function(x) x["samples"])
    cluster_samp <- cluster_samp[ind_leaf]
    cluster_samp <- lapply(cluster_samp, function(x) intersect(x, ovlap_subj))

    ind_not_zero <- which(unlist(sapply(cluster_samp,
                                        function(x) length(x))) != 0)
    cluster_samp <- cluster_samp[ind_not_zero]

    cluster_samp_unique <- lapply(cluster_samp, function(x) t(combn(x, 2)))
    df_1 <- distinct(as.data.frame(do.call("rbind", cluster_samp_unique)))

    is_leaf <- sapply(nodes2, function(x) x["is_leaf"])
    ind_leaf <- which(unlist(is_leaf) == TRUE)
    cluster_samp <- sapply(nodes2, function(x) x["samples"])
    cluster_samp <- cluster_samp[ind_leaf]
    cluster_samp <- lapply(cluster_samp, function(x) intersect(x, ovlap_subj))

    ind_not_zero <- which(unlist(sapply(cluster_samp,
                                        function(x) length(x))) != 0)
    cluster_samp <- cluster_samp[ind_not_zero]

    cluster_samp_unique <- lapply(cluster_samp, function(x) t(combn(x, 2)))
    df_2 <- distinct(as.data.frame(do.call("rbind", cluster_samp_unique)))

    df_3 <- data.frame(V1 = df_2$V2, V2 = df_2$V1)

    a <- intersect(df_1, df_2)
    b <- intersect(df_1, df_3)
    c <- rbind(a, b)

    numerator <- dim(df_1)[1] + dim(df_2)[1] - 2 * dim(c)[1]
    numerator <- numerator * obtain_d0_distance(tree1, tree2, data)
    denominator <- choose(length(ovlap_subj), 2)

    return(numerator/denominator)
}

#' Obtain d2 distance between two trees. (Deprecated)
#' @tree1 casual tree object 1.
#' @tree2 casual tree object 1.
#' @data the data used to train the casual forest.
obtain_d2_distance <- function(tree1, tree2, data) {
    d2 <- 0
    for (i in 1:length(rownames(data))) { # check all pair of data.
        y_hat1 <- predict_tree(tree1$nodes[[1]], data[i, ], tree1)[1]
        y_hat2 <- predict_tree(tree2$nodes[[1]], data[i, ], tree2)[1]
        d2 <- d2 + (y_hat1 - y_hat2)^2/length(rownames(data))
    }
    return(d2)
}

#' Obtain d2 distance between two trees (multicore edition).
#' @tree1 casual tree object 1.
#' @tree2 casual tree object 1.
#' @data the data used to train the casual forest.
obtain_d2_distance_mc <- function(tree1, tree2, data) {
    d2 <- foreach(i = 1:length(
        rownames(data)), .combine = "+") %dopar% { # check all pair of data.
        y_hat1 <- predict_tree(tree1$nodes[[1]], data[i, ], tree1)[1]
        y_hat2 <- predict_tree(tree2$nodes[[1]], data[i, ], tree2)[1]
        (y_hat1 - y_hat2)^2
    }
    d2 <- d2/length(rownames(data))
    return(d2)
}

#' Make prediction using single causal tree.
#' @node node object in causal tree object.
#' @sample single sample from the provided dataset (single patient).
#' @tree Causal tree object.
predict_tree <- function(node, sample, tree) {
  # Check if this is leaf node
  if (node$is_leaf == TRUE) {
    return(node$leaf_stats)
  } else {
    # Check split varibale
    split_var <- node$split_variable
    split_val <- node$split_value
    if (sample[split_var] <= split_val) {
      return(predict_tree(tree$nodes[[node$left_child]], sample, tree))
    } else {
      return(predict_tree(tree$nodes[[node$right_child]], sample, tree))
    }
  }
}

#' This is the function used for searching the representative tree.
#' @forestfile file name of the forest file, required type
#' "KIRC_ENSG00000162975_forest.RData".
#' @trainingset file name of the training data.
#' @n number of thread used.
#' @distance_type d0, d1, d1star, d2
get_reprtree <- function(forestfile, trainingset, n, distance_type = "d1") {
    # Load required packages.
    library(grf)
    library(randomForest)
    library(data.table)
    library(dplyr)
    library(doParallel)
    library(stringr)

    # Load a casual forest model
    load(forestfile) # This will create a variable "tau.forest" by default.

    # Load the training data.
    # Read the Training dataset.
    kirc_data <- data.table::fread(
        file = trainingset, header = TRUE) %>% as.data.frame()
    rownames(kirc_data) <- kirc_data$donorId
    kirc_data$donorId <- NULL

    # Obtain the W gene based on the given RData file.
    s <- forestfile
    gene_w <- stringr::str_match(s, "(?:_)(.*)(_)")[2] # Modify if needed.
    remove_w_index <- which(colnames(kirc_data) == gene_w)

    # Remove the W treatment effect column
    kirc_data <- kirc_data[, -c(remove_w_index)]

    # Get tree number.
    # tau.forest is from the loading tree object.
    tree_num <- tau.forest$`_num_trees`

    # Initialize the distance matrix
    d_matrix <- matrix(, nrow = tree_num, ncol = tree_num)

    # Calculate D(T).
    doParallel::registerDoParallel()
    options(cores = n)

    if (distance_type == "d2") {
        for (i in 1:tree_num) {
            mount_tree <- grf::get_tree(tau.forest, i)
            for (j in 1:tree_num) {
                if (i <= j) { # Same tree will be omitted.
                    next
                } else {
                    tmp_tree <- grf::get_tree(tau.forest, j)
                    dt_ij <- obtain_d2_distance_mc(mount_tree, 
                                                   tmp_tree, kirc_data)
                    d_matrix[i, j] <- dt_ij
                }
            }
        }
    } else if (distance_type == "d1") {
        for (i in 1:tree_num) {
            mount_tree <- grf::get_tree(tau.forest, i)
            for (j in 1:tree_num) {
                if (i <= j) { # Same tree will be omitted.
                    next
                } else {
                    tmp_tree <- grf::get_tree(tau.forest, j)
                    dt_ij <- obtain_d1_distance(mount_tree, tmp_tree)
                    d_matrix[i, j] <- dt_ij
                }
            }
        }
    } else if (distance_type == "d1star") {
        for (i in 1:tree_num) {
            mount_tree <- grf::get_tree(tau.forest, i)
            for (j in 1:tree_num) {
                if (i <= j) { # Same tree will be omitted.
                    next
                } else {
                    tmp_tree <- grf::get_tree(tau.forest, j)
                    dt_ij <- obtain_d1_star_distance(mount_tree,
                                                     tmp_tree, kirc_data)
                    d_matrix[i, j] <- dt_ij
                }
            }
        }
    } else if (distance_type == "d0") {
        for (i in 1:tree_num) {
            mount_tree <- grf::get_tree(tau.forest, i)
            for (j in 1:tree_num) {
                if (i <= j) { # Same tree will be omitted.
                    next
                } else {
                    tmp_tree <- grf::get_tree(tau.forest, j)
                    dt_ij <- obtain_d0_distance(mount_tree, tmp_tree, kirc_data)
                    d_matrix[i, j] <- dt_ij
                }
            }
        }
    }

    # fill the matrix
    for (i in 1:tree_num) {
        for (j in 1:tree_num) {
            if (i == j) {
                d_matrix[i, j] <- 0.0 # Distance between same tree should be 0.
            } else if (i > j) {
                next
            } else {
                d_matrix[i, j] <- d_matrix[j, i]
            }
        }
    }

    # Choose the smallest D(T) values.
    r <- colMeans(d_matrix)
    index_tree <- which.min(r)
    return(index_tree)

    #' Return the tree object. (Modify if needed)
    #' reptree <- get_tree(tau.forest, index_tree)
    #' return(reptree)
}