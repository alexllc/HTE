library(TCGAbiolinks)

# Fetch survival data from GDC
query <- GDCquery(project = paste0("TCGA-", project),
                data.category = "Clinical",
                file.type = "xml") # Obtain list from TCGAbiolinks:::getProjectSummary(<enter_project>)

#GDCdownload(query)
load(paste0("./clinical/", project, "_clinical.rda"))
# patient = GDCprepare_clinic(query, clinical.info = "patient")

#save(patient, file = paste0("./clinical/", project, "_clinical.rda"))
s.patient <-  c("bcr_patient_barcode", "gender","vital_status","days_to_birth", "days_to_death", "days_to_last_followup","race_list", "stage_event_pathologic_stage")

s.df <- NA
for (i in s.patient) {
    s.df <- cbind(s.df, patient[grep(i, colnames(patient))])
}
s.df <- s.df[,-1]
df_patient <- s.df[!duplicated(s.df),]

# Process patient info
df_patient$vital_status <- sapply(as.numeric(df_patient$vital_status), function(x) x - 1)
df_patient <- df_patient %>% mutate(survival_time = coalesce(days_to_death, days_to_last_followup)) # there are some negative days
max.censored <- max(df_patient$survival_time[df_patient$vital_status == 0])
df_patient <- df_patient[df_patient$survival_time >= 0,]
df_patient$vital_status[df_patient$survival_time == max.censored] <- 1
df_patient$imputed.log.times <- impute.survival(df_patient$survival_time, df_patient$vital_status)
df_patient$age <- floor(df_patient$days_to_birth/365)/-1# Convert age

extra <- c("days_to_birth", "survival_time", "days_to_last_followup", "days_to_death")
df_patient <- dplyr::select(df_patient, -extra)
df_patient <- df_patient[df_patient$imputed.log.times >= 0,]

df_patient[df_patient==""] <- NA

for (c in colnames(df_patient)) {

    if (!is.numeric(df_patient[,c]) && c != "bcr_patient_barcode") {
        which.one <- which( levels(df_patient[,c]) == "")
        levels(df_patient[,c])[which.one] <- NA
        df_patient[,c] = sapply(sapply(df_patient[,c], as.factor), as.numeric) 
        print(paste0(c, " is altered")) 
    }
}

df_patient = na.omit(df_patient)
