#Install RGLab packages
#devtools::install_github("DillonHammill/CytoExploreR")
#devtools::install_github("RGLab/flowCore",force=TRUE)
#devtools::install_github("RGLab/ggcyto",force=TRUE)
#devtools::install_github("RGLab/flowWorkspace",force=TRUE)
#devtools::install_github("RGLab/cytolib",force=TRUE)
#devtools::install_github("RGLab/CytoML",force=TRUE)
#devtools::install_github("RGLab/openCyto",force=TRUE)
#devtools::install_github("rlbarter/superheat")

#Load necessary packages
library(ggrepel)
library(flowCore)
library(tidyverse)
library(rsvd)
library(cytofWorkflow)
library(ComplexHeatmap)
library(Rtsne)
library(hash)
library(cowplot)
library(devtools)
library(flowViz)
library(kableExtra)
library(remotes)
library(CytoML)
library(flowAI)
library(flowWorkspace)
library(ggcyto)
library(CytoExploreR)
library(scales)

#Set seed
set.seed(1234)

#increase memory and use cache of larger files
knitr::opts_chunk$set(cache=TRUE, warning=FALSE, message=FALSE, cache.lazy=FALSE)

#Set WD to the project folder (not R project folder)
setwd("/home/dany/Documents/myvault/70-79_Travail/73_Recherche_emploi_Data_Scientist/roche_flowcytometrydataanalyst_202401/onsite_interview_preparation/skills_assessment/skill_assessment/")

#Importing data with CytoExploreR
data_folder_path <- '../data/FlowRepository_FR-FCM-Z4GN_files/'
fs <- cyto_load(data_folder_path)
conditions <- c(simplify(lapply(c(79:84),function(.i){paste("postSpn",as.character(.i),sep="_")})),
             simplify(lapply(c(1:6),function(.i){paste("preSpn",as.character(.i),sep="_")})),
             simplify(lapply(c(67:72),function(.i){paste("saline",as.character(.i),sep="_")})))

cyto_details(fs)$conditions <- conditions
new_channels <- simplify(base::lapply(cyto_channels(fs),function(.s) str_replace_all(.s,'FJComp.','')))
channels_to_plot <- cyto_channels(fs)[1:30]

#Plotting the arcsinh transform versus log transform
x <- seq(-1000,200000,10)
base::plot(x,asinh(x/6000),col="green")
lines(x,x)
y <- seq(0,200000,10)
points(x,log(x/6000),col="red")

# Check data to understand scaling
cyto_plot(fs[1:6],title=conditions,channels=c("CD19","7AAD"),legend=FALSE,density_modal=TRUE,group_by="conditions")

#Setup my own simple asinh transform and test plotting
##Beware of the confusion between flowcore AND cytoWorkflow
trans.obj <- function(x){asinh(x/6000)}
invtrans.obj <- function(x){sinh(x)*6000}
trans <- scales::trans_new("asinh_trans",trans.obj,invtrans.obj)
transList <- transformerList(c(cyto_fluor_channels(fs)),trans)
fs <- cyto_transform(fs,trans=transList,channels=cyto_fluor_channels(fs))

xtick<-c(-1e4,0,1e4,1e5,1e6)

cyto_plot(fs[1:6],title=conditions,channels=c("CD19","7AAD"),legend=FALSE,group_by="conditions",axes_trans=transList,
          xlim=c(-1e4,1e6),ylim=c(-1e4,1e6))





# #Get fcs files
# data_folder_path <- './data/FlowRepository_FR-FCM-Z4GN_files/'
# fcs_files <- list.files(data_folder_path)
# fs <- read.flowSet(files=fcs_files, path=data_folder_path, truncate_max_range = FALSE)
# 
# #Get colnames in flow cytometry data
# col_names <- BiocGenerics::colnames(fs)
# row_names <- BiocGenerics::rownames(fs)
# 
# #Put data into a dataframe for visualization
# mydata <- data.frame(fsApply(fs,flowCore::exprs),'sample_names'=rep(sampleNames(fs),fsApply(fs,nrow))) %>% 
#   tidyr::extract(sample_names,'sample','ivCD45_2negBorTcells_([a-zA-Z]{1,}_[0-9]{1,}).fcs') %>% 
#   separate(sample,into=c('condition','replica'),sep='_')
# 
# #Replace all '.' with '-'
# new_colnames <- colnames(mydata)
# #Discard 'A' for Area for all fluo column names
# new_colnames <- base::lapply(new_colnames,function(.s) if(grepl('FJComp.',.s)){return(substr(.s,1,nchar(.s)-2))}else{return(.s)})
# #Delete 'FJComp-' from all column names and change specific channel names
# new_colnames <- base::lapply(new_colnames,function(.s) str_replace_all(.s,'FJComp.',''))
# new_colnames <- base::lapply(new_colnames,function(.s) str_replace_all(.s,'_',''))
# new_colnames <- base::lapply(new_colnames,function(.s) str_replace_all(.s,'Alexa.Fluor.','AF'))
# new_colnames <- base::lapply(new_colnames,function(.s) str_replace_all(.s,'APC.Fire.750','APCFire750'))
# new_colnames <- base::lapply(new_colnames,function(.s) str_replace_all(.s,'PE.Cy5.5','PECy55'))
# new_colnames <- base::lapply(new_colnames,function(.s) str_replace_all(.s,'PE.Vio770','PEVio770'))
# new_colnames <- base::lapply(new_colnames,function(.s) str_replace_all(.s,'eFluor.','e'))
# new_colnames <- base::lapply(new_colnames,function(.s) str_replace_all(.s,'PE.e610','PEe610'))
# new_colnames <- base::lapply(new_colnames,function(.s) str_replace_all(.s,'PerCP.e710','PerCPe710'))
# new_colnames <- base::lapply(new_colnames,function(.s) str_replace_all(.s,'7AAD','SevenAAD'))
# #new_colnames <- base::lapply(new_colnames,function(.s) str_replace_all(.s,'.e','e'))
# 
# #Import the list of all relevant markers and probes
# marker_probe_df <- readr::read_csv("./ressources/markers_and_probes.csv") %>% 
#   mutate(new_colname=paste(marker,probe,sep=".")) %>% 
#   dplyr::rename('old_colname'=probe)
# 
# old <- marker_probe_df$old_colname
# new <- marker_probe_df$new_colname
# new_old_names_hash <- hash::hash(old,new)
# 
# # Finally change columns names
# new_colnames <- base::lapply(new_colnames,function(.s) if(.s %in% keys(new_old_names_hash)){return(new_old_names_hash[[.s]])}else{return(.s)})
# colnames(mydata) <- new_colnames
# 
# #Set up a proxy replica number to facilitate plotting
# replica_new <- as.character(rep(c(1:6),times=3))
# replica_old <- as.character(c(c(1:6),c(67:72),c(79:84)))
# 
# replica_df <- tibble(replica_old=replica_old,replica=replica_new)
# mydata <- mydata %>% dplyr::rename('replica_old'='replica') %>% 
#   left_join(replica_df,by=c('replica_old'))
# 
# 
# 
# 
# 
