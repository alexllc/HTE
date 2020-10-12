tr_ls = c()
for (i in 1:2000) {
ctree = get_tree(tau.forest, i)
if (length(ctree$nodes) > 3)
  tr_ls = c(i, tr_ls)
}

for (tree in tr_ls) {
  png(paste0("tree_", tree, "_plot.png"))
  plot(get_tree(tau.forest, tree))
  dev.off()
}