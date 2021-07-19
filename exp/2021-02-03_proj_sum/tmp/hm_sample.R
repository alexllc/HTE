data <- replicate(20, rnorm(50))
rownames(data) <- paste("Gene", c(1:nrow(data)))
colnames(data) <- paste("Sample", c(1:ncol(data)))


metadata <- data.frame(
      c(rep("case", ncol(data)/2), rep("control", ncol(data)/2)),
      c(rep("cond1", ncol(data)/4), rep("cond2", ncol(data)/4), rep("cond3", ncol(data)/4), rep("cond4", ncol(data)/4)),
      row.names=colnames(data))
colnames(metadata) <- c("casecontrol","condition")

metadata_gene <- data.frame(
      c(rep("Tcell", nrow(data)/2), rep("Bcell", nrow(data)/2)),
      row.names=rownames(data))
colnames(metadata_gene) <- c("Cell")


ann_colors = list(Time = c("white", "firebrick"),CellType = c(CT1 = "#1B9E77", CT2 = "#D95F02"),GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E"))