library(igraph)

GraphFeatures <- function(G_list,nrep){

  Global_summary <- NULL
  
  network_size=vector(mode = "list", length(G_list))
  FiedlerValue=vector(mode = "list", length(G_list))
  Adj_val=vector(mode = "list", length(G_list))
  Adj=vector(mode = "list", length(G_list))
  network_size=vector(mode = "list", length(G_list))
  most_infected_node=vector(mode = "list", length(G_list))
  mod=vector(mode = "list", length(G_list))
  deg=vector(mode = "list", length(G_list))
  cent=vector(mode = "list", length(G_list))
  trans=vector(mode = "list", length(G_list))
  di=vector(mode = "list", length(G_list))
  M=vector(mode = "list", length(G_list));spec=vector(mode = "list", length(G_list)); 
  Eigen_central=vector(mode = "list", length(G_list));
  most_infected_nod=vector(mode = "list", length(G_list))
  Rnot=vector(mode = "list", length(G_list));wc=vector(mode = "list", length(G_list));
  lc=vector(mode = "list", length(G_list));Adjacency_spectrum=vector(mode = "list", length(G))
  
for (g in 1:length(G_list)){
  for(rep in 1:nrep){   
    network_size[[g]][[rep]]<-vcount(G_list[[g]])
    #Adjacency matrix
    Adj[[g]][[rep]]<-as.matrix(G_list[[g]][])
    #Adjacency spectrum
    Adjacency_spectrum[[g]][[rep]]<-eigen(Adj[[g]][[rep]])
      
    
    M[[g]][[rep]]<- laplacian_matrix(G_list[[g]],sparse = FALSE)
    
    spec[[g]][[rep]]<- eigen(M[[g]][[rep]])

    FiedlerValue[[g]][[rep]]<- spec[[g]][[rep]]$values[length(spec[[g]][[rep]]$value)-1]
      
    
    Adj_val[[g]][[rep]]<-Adjacency_spectrum[[g]][[rep]]$values[1]
    
    Eigen_central[[g]][[rep]]<- Adjacency_spectrum[[g]][[rep]]$vectors[,ncol(Adjacency_spectrum[[g]][[rep]]$vectors)[1]]
  
    most_infected_node[[g]][[rep]]<- match(max(Eigen_central[[g]][[rep]]),Eigen_central[[g]][[rep]])# highest node degree

    Rnot[[g]][[rep]]<-(1/Adj_val[[g]][[rep]]) # spectral of adjacency for undirected network controls the epidemic threshold

    wc[[g]][[rep]]<- cluster_walktrap(G_list[[g]])  
    lc[[g]][[rep]]<- cluster_louvain(G_list[[g]]) 
    
    mod[[g]][[rep]]<- modularity(G_list[[g]], membership(lc[[g]][[rep]]))
     
    deg[[g]][[rep]]<- mean(degree(G_list[[g]], mode="in"))
      
   # cent[[g]][[rep]]<-  eigen_centrality(G_list[[g]], directed=T, weights=NA)
   # cent[[g]][[rep]] <- cent[[g]]$value
    
  #  trans[[g]]<- as.data.frame(transitivity(G_list[[g]], type="global"))
     
    di[[g]][[rep]]<- diameter(G_list[[g]], directed=F, weights=NA)
    
    
    Global_summary[[g]]<- as.data.frame(cbind(FiedlerValue[[g]], Adj_val[[g]], network_size[[g]], most_infected_node[[g]], deg[[g]], di[[g]])) 
  }
} 
    Global_summary <- do.call(rbind,  Global_summary)
#    Global_summary  <- as.data.frame(Global_summary)
  
    colnames(Global_summary) <- c('Fiedler','Adj_val','network_size','most_infected_node', 'Mean_degree','diameter')

   return(Global_summary)
}


