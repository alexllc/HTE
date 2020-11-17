# Using GRF to study HTE

## All project wide function scripts must be run from the base directory
`setwd("~/project/HTE/")`

To source these function scripts:

```
bin_ls = list.files("./bin")
for (bin in bin_ls){
    source(paste0("./bin/", bin))
}
```
