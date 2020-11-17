# This scripte is used to process the COVID data of UKBB from kenneth
# Created by Kai (Nov-15-2020)

require(data.table)

data <- load('/home/kenneth/out/UKBB/covid19/run_03Nov2020/without_imputation/df_cohorts_3Nov2020.RData')
s <- cbind(as.numeric(rownames(df.a)), df.a[,c(1,2,4:6,8:9,11,16:21,28,30,34,36,76:80,82,90,97,98)])
colnames(s)[c(1,28)] <- c('eid', 'Severity')
rownames(s) <- NULL


d <- cbind(as.numeric(rownames(df.b)), df.b[,c(1,2,4:6,8:9,11,16:21,28,30,34,36,76:80,82,90,97,98)])
colnames(d)[c(1,28)] <- c('eid', 'U071')
rownames(d) <- NULL

col_names <- colnames(s) 
ord_names <- c(col_names[1], col_names[c(-1,-28)][order(col_names[c(-1, -28)])], col_names[28]) 
s <- s[,ord_names]
colnames(s)[c(-1, -28)] <-
c('AF','21022-0.0','20117-0.1','22127-0.0','21001-0.0','30710-0.0','CAD','22130-0.0','30700-0.0','Depression','4079-0.0','2453-0.0','21000-0.1','30750-0.0','30760-0.0','Heart_failure','Hypertension','30780-0.0','31-0.0','20116-0.1','4080-0.0','2443-0.0','2443-1.0','30870-0.0','30670-0.0','whr')

col_names <- colnames(d) 
ord_names <- c(col_names[1], col_names[c(-1,-28)][order(col_names[c(-1, -28)])], col_names[28]) 
d <- d[,ord_names]
colnames(d)[c(-1, -28)] <-
c('AF','21022-0.0','20117-0.1','22127-0.0','21001-0.0','30710-0.0','CAD','22130-0.0','30700-0.0','Depression','4079-0.0','2453-0.0','21000-0.1','30750-0.0','30760-0.0','Heart_failure','Hypertension','30780-0.0','31-0.0','20116-0.1','4080-0.0','2443-0.0','2443-1.0','30870-0.0','30670-0.0','whr')

write.csv(s, '~/data/UKBB/UKBB_clinical_variable_severity_with_imputation.csv', row.names = F)
write.csv(d, '~/data/UKBB/UKBB_clinical_variable_death_with_imputation.csv', row.names = F)



colnames(data)[c(4,14,21)] <- c('20117-0.0','21000-0.0','20116-0.0')

idx <- which(colnames(data) == '20117-0.0')
c1 <- as.numeric(data[,idx, with = F] == 0)
c2 <- as.numeric(data[,idx, with = F] == 1)
c3 <- as.numeric(data[,idx, with = F] == 2)

cls <- cbind(c1, c2, c3)
colnames(cls) <- c('20117-0.1', '20117-0.2', '20117-0.3')
cls[is.na(data[,idx, with = F]),] <- NA

data <- cbind(data[,1:(idx-1), with = F], cls, data[,(idx+1):ncol(data), with = F])



idx <- which(colnames(data) == '20116-0.0')
c1 <- as.numeric(data[,idx, with = F] == 0)
c2 <- as.numeric(data[,idx, with = F] == 1)
c3 <- as.numeric(data[,idx, with = F] == 2)

cls <- cbind(c1, c2, c3)
colnames(cls) <- c('20116-0.1', '20116-0.2', '20116-0.3')
cls[is.na(data[,idx, with = F]),] <- NA
data <- cbind(data[,1:(idx-1), with = F], cls, data[,(idx+1):ncol(data), with = F])



idx <- which(colnames(data) == '21000-0.0')
c1 <- as.numeric(data[,idx, with = F] == 1)
c2 <- as.numeric(data[,idx, with = F] == 2)
c3 <- as.numeric(data[,idx, with = F] == 3)
c4 <- as.numeric(data[,idx, with = F] == 4)
c5 <- as.numeric(data[,idx, with = F] == 5)
c6 <- as.numeric(data[,idx, with = F] == 6)

cls <- cbind(c1, c2, c3, c4, c5, c6)
colnames(cls) <- c('21000-0.1', '21000-0.2', '21000-0.3','21000-0.4', '21000-0.5', '21000-0.6')
cls[is.na(data[,idx, with = F]),] <- NA

data <- cbind(data[,1:(idx-1), with = F], cls, data[,(idx+1):ncol(data), with = F])

