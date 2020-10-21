library(igraph)
library(EpiModel)
library(intergraph)
library(DiagrammeR)
library(tidyverse)

source( "epi_sim.R")

source( "epi-summary.R")

source( "network_summary.R")


Network_Features <- function(models,numb_simult){
  

  network_summary <- NULL
  
  g=list()
  network_size=list()
  g_names=list()
  avg_deg=list()
  
  for (i in 1:length(models))
  {
    g[[i]] <- models[[i]]
    
    g_names[[i]] <- names(models[i])
  
  
    network_size[[i]]<-vcount(g[[i]])
    
    avg_deg[[i]] <- mean(degree(g[[i]], mode="all", normalized=TRUE))


    network_summary[[i]] <- c(g_names[[i]], network_size[[i]], avg_deg[[i]]) 
    
  }
  
  
  network_summary <- do.call(rbind, network_summary)
  network_summary  <- as.data.frame(network_summary)
  
  colnames(network_summary) <- c('Network_type', 'Network_size','Average_degree')

  
  
  return(network_summary)
}

# g
# names(g) <- c("random_graph","small_world","lattice","scale_free","spatial")
# numb_simult=5
# g1= NULL
# for (j in 1:numb_simult){
#   if(names(g)[[i]]=="random_graph"){
#     g1[[j]]=erdos(4,.5)}
#   
#   if(names(g)[[i]]=="small_world"){
#     g1[[j]]=sm_world(1, 4, 2, 0.05)}
# }
