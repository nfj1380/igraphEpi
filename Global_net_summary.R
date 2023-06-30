Global_net_summary <- function(folder_name){
  
  networks <- list.files(paste(folder_name))
  
  #}  
  #for (i in 1:length(networks))##
  Global_summary  <- foreach( i = 1:length(networks), combine=rbind,  .packages=c('tidyverse', 'igraph','intergraph', 'DiagrammeR')) %dopar% 
    
    {
      source('simPathE.R')
      source('countStates.R ') 
      
      print(networks[i])
      filename=paste(folder_name, networks[i], sep="/") #
      
      dat <- read.table(filename) 
      g <- igraph::graph_from_data_frame(dat) %>%  as.undirected( "collapse")
      
      
      #g <- read_graph(networks[4], format = c("edgelist")) %>% as.undirected( "collapse") #seem only to work on smaller graphs?
      #one edge for each pair of connect vertices
      g_names <- gsub(".edges","",networks[i]) ##
      
      netclust <- clusters(g) #look for subgraphs
      gcc <- V(g)[netclust$membership == which.max(netclust$csize)]#select vertices from the largest sub-graph
      #we could do this for other component graphs?
      gccnet <- induced.subgraph(g, gcc) #make it a igraph object.
      
      #need to remove loops etc
      gcc_s <- igraph::simplify(gccnet, remove.multiple = TRUE, remove.loops = TRUE)
      #plot(gcc_s,main=networks[1])##

#network characteristics
      
      #network_size - number of nodes
network_size<-vcount(gcc_s)
#convert to a laplacian matrix
    
      
M <- laplacian_matrix(gcc_s,sparse = FALSE)

#Adjacency matrix
Adj <-as.matrix(gcc_s[])
#Adjacency spectrum
Adjacency_spectrum<-eigen(Adj)

#calc eigenvectors
spec <- eigen(M)
#take the second smallest eigenvalue
FiedlerValue<- spec$values[length(spec$value)-1]

#Largest value of the adjacency
##This value controls the propagation of infection in SIS/SIR model
### The larger Adj_val the faster the growth of infection at the beginning time

Spectral_radius<-Adjacency_spectrum$values[1]

Eigen_central <- Adjacency_spectrum$vectors[,ncol(Adjacency_spectrum$vectors)[1]]

highest_degree <- match(max(Eigen_central),Eigen_central)# highest node degree

highest_degree_scaled <- highest_degree/network_size 

#Note:Probability of infection of nodes depends on V_1(first eigen vector of adjacency matrix)
## so the prob that a node is infected is proportional to its eigen vector centrality (V_1)
## since centrality depends on number of nodes its connected to and how well node are connected

#Newman's relative Q (from Sah et al 2017) modularity

#detect communities and calc relative Newman's Q
# wc <- cluster_walktrap(gcc_s) # This function tries to find densely connected subgraphs, 

lc <- cluster_louvain(gcc_s) #this is the method used inSah et al 2017 - allows comparison? 
#also called communities in a graph via random walks.

mod <- modularity(gcc_s, membership(lc))

max_mod <- assortnet::assortment.discrete(Adj, types = as.factor(lc$membership), weighted = T)$r 

#max_mod <- 
Qrel <- mod/max_mod #Sah et al normalize modularity this way, but Im not sure why maximum modulairity
#is the best way to do this

#other network characteristics
deg <- mean(degree(gcc_s, mode="all", normalized=TRUE)) #normalized currently

trans <- as.data.frame(transitivity(gcc_s, type="global")) #clustering coefficient

di <- diameter(gcc_s, directed=F, weights=NA)

mean_dist <- mean_distance(gcc_s, directed=F)


data.frame(Network=g_names, 
           Fiedler=  FiedlerValue, Spectral_radius=Spectral_radius, Network_size=network_size, 
           Highest_degree=highest_degree, highest_degree_scaled= highest_degree_scaled,
           Modularity=mod, Qrel =Qrel, Mean_degree=deg,
           Transivity=trans, Diameter=di, Mean_pathLength= mean_dist)
    }

  Global_summary <- do.call(rbind, Global_summary)
  
  return(Global_summary)
  
  parallel::stopCluster(cl)
}
  