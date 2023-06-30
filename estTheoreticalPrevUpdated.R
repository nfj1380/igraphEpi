estTheoreticalPrev<- function(folder_name, R0=R0){
  
  networks <- list.files(paste(folder_name))
  
  #Global_summary <- NULL
  
  #make sure all network object in the folder can load.
  
  netcheck <- function(networks ){ 
    filename=paste(folder_name, networks, sep="/") #
    
    dat <- read.table(filename) 
    options(show.error.messages = TRUE)
    print(filename)
    gtry <- try(igraph::graph_from_data_frame(dat) %>%  as.undirected( "collapse"))
    if(is( gtry,"try-error")) {print('graph error')}
    else { print('graph ok')}
  }
  
  # removing mammalia-voles-kcs-trapping-01 as there is a strange error with this one
  
  #will throw an error if there is anything wrong with any graph in the folder
  #}  
  #for (i in 1:length(networks))##
  Global_summary  <- foreach( i = 1:length(networks), combine=rbind,  .packages=c('tidyverse', 'igraph','intergraph', 'DiagrammeR')) %dopar% 
    
    {
      print(networks[i])#
      filename=paste(folder_name, networks[i], sep="/") #
      
      dat <- read.table(filename) #haven't dealt with the temporal nature of this yet. Third column appears to be a date...
      g <- igraph::graph_from_data_frame(dat) %>%  as.undirected( "collapse")
    
      g_names <- gsub(".edges","",networks[i])  #
      
    netclust <- clusters(g)
    gcc <- V(g)[netclust$membership == which.max(netclust$csize)]
    gccnet <- induced.subgraph(g, gcc)
    gcc_s <- igraph::simplify(gccnet, remove.multiple = TRUE, remove.loops = TRUE)
    
    degree_seq <- igraph::degree(gcc_s)
    
    max_degree <- max(degree_seq)
    p_k <- numeric(max_degree)
    
    for (i in 1:max_degree) {
      p_k[i] <- sum(degree_seq == i) / vcount(g)
    }
    
    S_0 <- 1.0 #ll vertices susceptible
    
    k <- mean(p_k)
    
    theta <- k/var(p_k)-k #to account for the inability of node to infect its parent node. Newman 2002

    # Compute the expected number of steps T_k required for an infection to reach a node of degree k
    
    S_inf_49 <- exp(-R0 * (1 - S_0)) 
    # S_inf_50 <- sum(p_k * (1 - exp(-R0 * k * T_k))) # includes degree structure
    S_inf_50 <- sum(p_k * (1 - exp(-R0 * k * theta)))
    
  #  return(list(S_inf_49 = S_inf_49, S_inf_50 = S_inf_50))
    
    data.frame(Network=g_names, 
               S_inf_49 = S_inf_49 ,
               S_inf_50  =  S_inf_50)
    }
  
  #write.table(file="globalData.csv",append=TRUE,sep=",",col.names=TRUE,row.names=TRUE)
  #Global_summary[[i]] 
  
  #Global_summary_list 
  #write.table(Global_summary[,,i],file="globalData.i.csv",append=TRUE,sep=",",col.names=TRUE,row.names=TRUE) 
  
  Global_summary <- do.call(rbind, Global_summary)
  
  
  # colnames(Global_summary) <- c('Network','Avg_infected_R02','Avg_infected_R05', 'Fiedler','Adj_val','Rnot','network_size','most_infected_node', 'Modularity', 'Mean_degree', 'Centrality', 'Transivity', 'Diameter')
  # Global_summary  <- as.data.frame(Global_summary)
  
  return(Global_summary)
  
  parallel::stopCluster(cl)
  
}