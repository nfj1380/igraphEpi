## Disease simulation


library(igraph)
#Simple SIR disease simulation
# Disease status: 0=healthy/susceptible, 1=infected, 2=removed/dead/recovered

episim <- function(Graphs, nticks=10, beta=0.1, gamma=0.2, propInfected){
  
  
# if (is.igraph(Graphs)){
#    Graphs=as.list(Graphs)
#    }
for (i in 1:length(Graphs)){
  if (is.igraph(Graphs[[i]])){
    Graphs[[i]]=as.list(Graphs[[i]])
    }
#    return(Graphs)
  }
  
  E=list()
  n=list()
  inf_per_time=list()
  infected=list()
  clear=list()
  stat=list()
  nextStat=list()
  
  
for (G in 1:length(Graphs)){
  E[[G]] <- as.matrix(get.adjacency(Graphs[[G]]))
  n[[G]] <- ncol(E[[G]]) #number of columns (nodes) in the list of network object
  inf_per_time[[G]]<-matrix(data=0,nticks+1,n[[G]])#list of matrix related to the igraph objects with progression of infection/transition of state of each node for each time step
  
    
#  }
    
#  result=list()
 # graph_idx=1
  
#  for (G in Graphs){
    ## G is an igraph object
#    E <- get.adjacency(G,type = c("both")) #adjacency matrix of graph (network object)
 #   n <- ncol(E) #number of columns (nodes) in network object
#    inf_per_time<-matrix(data=0,nticks+1,n)#progression of infection/transition of state of each node for each time step
    
    # Set initial state  
    # initState is list of initial state of each vertex in G.  E.g., from (1,0,0...0) or (0,1,0,1,0,....)
    
    # if (is.null(initState[[G]])) {
    #   infected[[G]] <- rep(1,numInfected)
    #   clear[[G]]<-rep(0,(n[[G]]-numInfected))
    #   stat[[G]]<-sample(c(infected[[G]], clear[[G]]), n[[G]], replace=FALSE)
    #   
      
    # } else {
    #   stat[[G]] <- initState[[G]] #initial existing state.
    # }
    # 
  
     stat[[G]] <- rbinom(n[[G]],1,propInfected)#set initial state of each nodes if non (note, no recovered node at initial state)
    
     current_stat=stat
    
     nextStat[[G]]<-rep(0,length(stat[[G]]))
     time=2 # since we have nticks+1
    
    #note: there is no recovered node at initial state only susceptible nodes and some infected nodes
    #initial state vector must not contain 2, but only 0 and 1 for the first row
    #set next state vector initially as initial state.
    for (tick in 1:nticks) {
      inf_per_time[[G]][1,]=current_stat[[G]]#time 0=initial state
      # first transition phase: S -> I
      for (i in which(stat[[G]]==0)) {#for all susceptible nodes (denoted 0) in initial state
        for (j in (1:n[[G]])) { # for all other neighbor nodes
          if ((E[[G]][i,j]>=1)&(stat[[G]][j]==1)&(runif(1)<=beta)) #each infected node infects its nearest neighbor/node it shares a link/edge with, with prob beta
          {
            nextStat[[G]][i] <- 1;break#assign node as infected if above condition is met else and break out of loop
          }
          else{
            nextStat[[G]][i]<-0 #node remains as susceptible, since not all contact leads to an infection.
          }
        }
      }
      
      # second transition phase: I -> R
      for (i in which(stat[[G]]==1)) {#for all infected nodes in initial state
        if (runif(1)<=gamma) { #compares a randomly generated uniform number to recovery rate
          nextStat[[G]][i] <- 2 #assigns node as recovered if condition above is met else
        }
        else{
          nextStat[[G]][i]<-1 #nodes remain infected 
        }
      }
      
      #for SIR, recovered node do not participate in further disease propagation (permanent immunity)    
      
      stat[[G]] <- nextStat[[G]]#assigns states of nodes for each time step as existing state 
      
      inf_per_time[[G]][time,] <- stat[[G]]
      rownames(inf_per_time[[G]])=0:nticks#re-order row names from time step 0=initial state
      colnames(inf_per_time[[G]])=NULL#make column names as null
      
      
        time=time+1
      
     }
    
   } 
  return(inf_per_time)
}



# G=f1
# for (i in 1:length(G)){
#   if (is.igraph(G[[i]])){
#     G[[i]]=as.list(G[[i]])
#   }
#   return(G)
# }
# 
# E=list()
# n=list()
# for (k in (1:length(G))){
# E[[k]] <- as.matrix(get.adjacency(G[[k]]))
# n[[k]] <- ncol(E[[k]]) #number of columns (nodes) in network object
# }
