#Packages to install
library(ggrepel)
library(flowCore)
library(tidyverse)
library(rsvd)
#library(cytofkit2) deprecated
library(cytofWorkflow)
library(ComplexHeatmap)
library(Rtsne)
library(hash)

#Set seed
set.seed(1234)

#Set WD to the project folder (not R project folder)
setwd("/home/dany/Documents/myvault/70-79_Travail/73_Recherche_emploi_Data_Scientist/roche_flowcytometrydataanalyst_202401/onsite_interview_preparation/skills_assessment/")

#Get fcs files
data_folder_path <- './data/FlowRepository_FR-FCM-Z4GN_files/'
fcs_files <- list.files(data_folder_path)
fs <- read.flowSet(files=fcs_files, path=data_folder_path, truncate_max_range = FALSE)

#Get colnames in flow cytometry data
col_names <- BiocGenerics::colnames(fs)
row_names <- BiocGenerics::rownames(fs)

#Put data into a dataframe for visualization
mydata <- data.frame(fsApply(fs,flowCore::exprs),'sample_names'=rep(sampleNames(fs),fsApply(fs,nrow))) %>% 
  tidyr::extract(sample_names,'sample','ivCD45_2negBorTcells_([a-zA-Z]{1,}_[0-9]{1,}).fcs') %>% 
  separate(sample,into=c('condition','replica'),sep='_')

#Replace all '.' with '-'
new_colnames <- colnames(mydata)
#Discard 'A' for Area for all fluo column names
new_colnames <- base::lapply(new_colnames,function(.s) if(grepl('FJComp.',.s)){return(substr(.s,1,nchar(.s)-2))}else{return(.s)})
#Delete 'FJComp-' from all column names and change specific channel names
new_colnames <- base::lapply(new_colnames,function(.s) str_replace_all(.s,'FJComp.',''))
new_colnames <- base::lapply(new_colnames,function(.s) str_replace_all(.s,'_',''))
new_colnames <- base::lapply(new_colnames,function(.s) str_replace_all(.s,'Alexa.Fluor.','AF'))
new_colnames <- base::lapply(new_colnames,function(.s) str_replace_all(.s,'APC.Fire.750','APCFire750'))
new_colnames <- base::lapply(new_colnames,function(.s) str_replace_all(.s,'PE.Cy5.5','PECy55'))
new_colnames <- base::lapply(new_colnames,function(.s) str_replace_all(.s,'PE.Vio770','PEVio770'))
new_colnames <- base::lapply(new_colnames,function(.s) str_replace_all(.s,'eFluor.','e'))
new_colnames <- base::lapply(new_colnames,function(.s) str_replace_all(.s,'PE.e610','PEe610'))
new_colnames <- base::lapply(new_colnames,function(.s) str_replace_all(.s,'PerCP.e710','PerCPe710'))
new_colnames <- base::lapply(new_colnames,function(.s) str_replace_all(.s,'7AAD','SevenAAD'))
#new_colnames <- base::lapply(new_colnames,function(.s) str_replace_all(.s,'.e','e'))

#Import the list of all relevant markers and probes
marker_probe_df <- readr::read_csv("./ressources/markers_and_probes.csv") %>% 
  mutate(new_colname=paste(marker,probe,sep=".")) %>% 
  dplyr::rename('old_colname'=probe)

old <- marker_probe_df$old_colname
new <- marker_probe_df$new_colname
new_old_names_hash <- hash::hash(old,new)

# Finally change columns names
new_colnames <- base::lapply(new_colnames,function(.s) if(.s %in% keys(new_old_names_hash)){return(new_old_names_hash[[.s]])}else{return(.s)})
colnames(mydata) <- new_colnames

#Set up a proxy replica number to facilitate plotting
replica_new <- as.character(rep(c(1:6),times=3))
replica_old <- as.character(c(c(1:6),c(67:72),c(79:84)))

replica_df <- tibble(replica_old=replica_old,replica=replica_new)
mydata <- mydata %>% dplyr::rename('replica_old'='replica') %>% 
  left_join(replica_df,by=c('replica_old'))





