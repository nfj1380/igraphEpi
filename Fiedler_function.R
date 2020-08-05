
#-------------------------------------------
# Calculating network properties for multiple graphs
#-------------------------------------------
#code by Nick Fountain-Jones May/June 2020

library(igraph)
library(EpiModel)
library(intergraph)
library(DiagrammeR)
library(tidyverse)


#this function calculates prevalence ofan SIR pathogen network as well as the corresponding network characteristics
Network_sum <- function(folder_name){
  
  networks <- list.files(paste(folder_name))
  
  Global_summary <- NULL
  
  
  for (i in 1:length(networks))
  {
    print(networks[i])
    filename=paste(folder_name, networks[i], sep="/")
    dat <- read.table(filename) #haven't dealt with the temporal nature of this yet. Third column appears to be a date...
    g <- graph_from_data_frame(dat)%>% as.undirected( "collapse")
    
    #g <- read_graph(networks[4], format = c("edgelist")) %>% as.undirected( "collapse") #seem only to work on smaller graphs?
    #one edge for each pair of connect vertices
    g_names <- gsub(".edges","",networks[i])
    
    netclust <- clusters(g) #look for subgraphs
    gcc <- V(g)[netclust$membership == which.max(netclust$csize)]#select vertices from the largest sub-graph
    #we could do this for other component graphs?
    gccnet <- induced.subgraph(g, gcc) #make it a igraph object.
    
    #need to remove loops etc
    gcc_s <- igraph::simplify(gccnet, remove.multiple = TRUE, remove.loops = TRUE)
    plot(gcc_s,main=networks[i])
    
    #------------------------------------------
    #SIR model from igraph. 
    
    sm2 <- sir(gcc_s, beta=2, gamma=1, no.sim=1000)#R0 of 2 #crashes when I do decimals for beta/gamma
    sm5 <- sir(gcc_s, beta=5, gamma=1, no.sim=1000)#R0 of 5. Not sure how the number of timesteps are calculated
    
    time_bins(sm2, middle = TRUE) #sets the time bins to calculate stats
    time_bins(sm5, middle = TRUE)
    
    med2 <- as.data.frame(median(sm2))
    med2$NI[is.na(med2$NI)] <- 0 #not sure about why Im getting NAs here (slightly concerned - but making them 0s ata this stage) 
    med5 <- median(sm5)
    med5$NI[is.na(med5$NI)] <- 0
    #make sure NAs are removed
    #plot(sm2)
    #plot(sm5)
    propI2 <- as.data.frame((mean(med2$NI))/length(V(gcc_s)))
    propI5 <- as.data.frame((mean(med5$NI))/length(V(gcc_s))) ##divides number infected over time by the number of nodes
    
    #we could make this maxium prevalence or any other characteristics. Need to think more about this.
    
    #------------------------------------------
    
    #network characteristics
    
    #convert to a laplacian matrix
    M <- laplacian_matrix(gcc_s,sparse = FALSE)
    
    #calc eigenvectors
    spec <- eigen(M)
    #take the second smallest eigenvalue
    FiedlerValue<- spec$values[length(spec$value)-1]
    
    #Newman's relative Q (from Sah et al 2017) modularity
    
    #detect commnuities and calc relative Newman's Q
    wc <- cluster_walktrap(gcc_s) # This function tries to find densely connected subgraphs, 
    lc <- cluster_louvain(gcc_s) #this is the method used inSah et al 2017 - allows comparison? 
    #also called communities in a graph via random walks.
    
    mod <- modularity(gcc_s, membership(lc))
    #max_mod <- 
    #Qrel <- mod/max_mod #Sah et al normalize modularity this way, but Im not sure why maximum modulairity is the best way to do this
    #plot(lc,gccnet)  check out the plot if needed
    
    #other network characteristics
    deg <- mean(degree(gcc_s, mode="all", normalized=TRUE)) #normalized currently
    
    cent <-  eigen_centrality(gcc_s, directed=F, weights=NA)
    cent <- cent$value
    
    trans <- as.data.frame(transitivity(gcc_s, type="global"))
    
    di <- diameter(gcc_s, directed=F, weights=NA)
    
    #Global_summary[1] <- bind_cols(g_names, propI2,propI5, FiedlerValue, mod, deg, cent, trans, di) %>%
    #set_names('g_names')
    
    Global_summary[[i]] <- c(g_names, propI2,propI5, FiedlerValue, mod, deg, cent, trans, di) 
  }
  
  Global_summary <- do.call(rbind, Global_summary)
  colnames(Global_summary) <- c('Network','Avg_infected_R02','Avg_infected_R05', 'Fiedler', 'Modularity', 'Mean_degree', 'Centrality', 'Transivity', 'Diameter')
  Global_summary  <- as.data.frame(Global_summary)
  
  
  return(Global_summary)
}


#create graph object
folder_name = c('Networks') #.edges file need to be in the working directory too at this stage.

test <- Network_sum(folder_name)

