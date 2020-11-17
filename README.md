# Heterogeneous treatment effect

## All binaries must be run from the base directory
`setwd("~/project/HTE/")`

## To source bin scripts
```
bin_ls = list.files("./bin")
for (bin in bin_ls){
    source(paste0("./bin/", bin))
}
```
