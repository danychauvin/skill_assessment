#Packages to install
library(ggrepel)
library(flowCore)
library(tidyverse)
library(rsvd)
#library(cytofkit2) deprecated
library(cytofWorkflow)
library(ComplexHeatmap)
library(Rtsne)

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

#Plot density histogram of FSC height
mydata %>% 
  select(FSC.H,condition,replica) %>%
  ggplot(aes(FSC.H))+
  geom_density(aes(col=replica))+
  facet_wrap(~condition)

#How many cells per condition
mydata %>% dplyr::count(condition,replica)



