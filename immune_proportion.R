# Before entering R session
# wget http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/LUSC/20160128/gdac.broadinstitute.org_LUSC.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0.tar.gz

library(EPIC)
library(data.table)
lusc = fread('LUSC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt')
dd = as.data.frame(lusc) # don't leave out!
dd = dd[,dd[1,]=="scaled_estimate"]
dd = as.matrix(dd)
mode(dd)='numeric'
dd=dd*1e6 
tmp=strsplit(lusc$`Hybridization REF`,'\\|')
as.character(tmp)
tmp=sapply(tmp,function(x)x[[1]])  ##extract the first value of the matrix
tmp.vv=which(nchar(tmp)>1)
rownames(dd) = tmp
dd = dd[tmp.vv,]
bulk <- dd[-1,] #remove the gene_id row
bulk <- bulk[!duplicated(rownames(bulk)),]
out = EPIC(bulk)
cells = out$cellFractions

# Only select normal tumor samples
# you can use these two lines to check again
tmp = strsplit(rownames(cells), '-')
table(sapply(tmp,function(x)x[[4]]))

cells = cells[grep("-01A-|-01B-", rownames(cells)),]
tmp = strsplit(rownames(cells), '-')
tmp = sapply(tmp, function(x) paste(x[[1]],x[[2]],x[[3]], sep='-'))
rownames(cells) = tmp
cells = as.data.frame(cells)
write.csv(cells, "LUSC_immune_cells.csv", row.names=T)