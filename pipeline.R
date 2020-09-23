#' Vignette to run the EpicGraph package

#Load packages - keep this to a minimum
library(igraph)
library(EpiModel)
library(intergraph)
library(DiagrammeR)
library(tidyverse)

# load all function codes. This will disappear when we formally make this a function
#source("getGraphFeatures.R")


#Nicks experiment

#-----------------------------------------------------------------------------------
#Step 1: Sumamrise network features
#-----------------------------------------------------------------------------------
#folder_name = c('Networks') #.edges file need to be in the working directory too at this stage.

#NetSum <- getGraphFeatures(network_name = f,  multinetwork = TRUE)

#NetSum

#Raima's code

source(knitr::purl("Experiment-1.Rmd", quiet=TRUE))

G <- G_2 # Sample graph to be given

# for (i in 1:length(G)){
#  if (is.connected(G[[i]])){
#    G=G
#  }else{
#    print("G is not connected")
#  }
# }

netsim<-episim(G,beta = 0.3,gamma=0.2) # epidemic simulation

netsummary<-epi_summary(netsim) # epidemic summary on a network

All_netsummary<-nets_epi_summary(G) # epidemic summary on all networks

All_netsummary

system.time({All_netsummary})

