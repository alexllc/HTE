library(readxl)

file_ls = list.files(".")

for (file in file_ls) {
    ifile = read_excel(file)
    fname = strsplit(file, "\\.")[[1]][1]
    write.csv(ifile, file = paste0(fname, ".csv"), row.names = FALSE)
}