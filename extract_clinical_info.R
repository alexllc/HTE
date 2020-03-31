library(TCGAbiolinks)

cancerType = "HNSC"
## Retreive clinical information (not survival) from TCGAbiolinks

query <- GDCquery(project = paste0("TCGA-", cancerType),
       				  data.category = "Clinical",
       				  file.type = "xml") # Obtain list from TCGAbiolinks:::getProjectSummary(<enter_project>)

		GDCdownload(query)

		clinical_set <- c("drug", "patient")

		for (c in clinical_set) {
		assign(c, GDCprepare_clinic(query, clinical.info = c))
		}