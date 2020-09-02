#' Vignette to run the EpicGraph package

#Load packages - keep this to a minimum
library(igraph)
library(EpiModel)
library(intergraph)
library(DiagrammeR)
library(tidyverse)

# load all function codes. This will disappear when we formally make this a function
source(".getGraphFeatures.R")


#-----------------------------------------------------------------------------------
#Step 1: Sumamrise network features
#-----------------------------------------------------------------------------------
folder_name = c('Networks') #.edges file need to be in the working directory too at this stage.

NetSum <- getGraphFeatures(network_name = folder_name,  multinetwork = TRUE)

NetSum