library(dplyr)
library(data.table)

donor = fread("donor.tsv")

col_selected = c("icgc_donor_id","project_code", "donor_sex",  "donor_vital_status",  "donor_age_at_diagnosis" , "donor_survival_time","donor_interval_of_last_followup" , "disease_status_last_followup")

clinical = filter(donor, !is.na(donor_survival_time) | !is.na(donor_interval_of_last_followup)) 

%>% select(all_of(col_selected))

mutation = fread("simple_somatic_mutation.open.tsv")