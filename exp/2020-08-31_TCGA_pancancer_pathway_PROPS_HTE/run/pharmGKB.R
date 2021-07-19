library("rBiopaxParser")

# PharmGKB
# Downloaded BioPax XML from https://www.pharmgkb.org/downloads

pharmgkb_loc <- "~/project/HTE/raw/pharmgkb/"
pharmgkb_files <- list.files(pharmgkb_loc)

bpax <- pharmgkb_files[grep("*.owl", pharmgkb_files)]

for (bpx_file in bpax) {
    biopax <- readBiopax(paste0(pharmgkb_loc, bpx_file))
    print(biopax)
    pw_list = listInstances(biopax, class="pathway")
    pw_complete = selectInstances(biopax, class="pathway")
    pwid1 = ""
    getInstanceProperty(biopax, pwid1, property="NAME")
    pw_1_component_list = listPathwayComponents(biopax,pwid1)
    pw_1_components = selectInstances(biopax,id=pw_1_component_list$id)
    pw_1_adj = pathway2AdjacancyMatrix(biopax, pwid1, expandSubpathways=TRUE, splitComplexMolecules=TRUE, verbose=TRUE)
}