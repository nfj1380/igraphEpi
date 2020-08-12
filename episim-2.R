##

getGraphFeatures <- function() {
  ## output: data frame
  
}

# converting downloaded networks into igraph objects
# use something liek this to convert:
#dat <- read.table(networks[i]) #haven't dealt with the temporal nature of this yet. Third column appears to be a date...
#g <- graph_from_data_frame(dat)%>% as.undirected( "collapse")

G <- graph_from_data_frame(dat)%>% as.undirected( "collapse")

library(igraph)
#Simple SIR disease simulation
# Disease status: 0=healthy/susceptible, 1=infected, 2=removed/dead/recovered
episim <- function(G, nticks=10, beta=0.1, gamma=0.2, initState=NULL, propInfected=0.5) {
  ## G is igraph object
  E <- G[] #adjacency matrix of graph (network object)
  n <- ncol(E) #number of columns (nodes) in network object
  inf_per_time<-data.frame(matrix(0,nticks+1,n))#progression of infection/transition of state of each node for each time step

  # Set initial state  
  # initState is list of initial state of each vertex in G.  E.g., from (1,0,0...0) or (0,1,0,1,0,....)
  if (is.null(initState)) {
    stat <- rbinom(n,1,propInfected)#set initial state of each nodes if non (note, no recovered node at initial state)
  } else {
    stat <- initState #initial existing state.
  }
  current_stat=stat
  nextStat<-rep(0,length(stat))
  time=2 # since we have nticks+1
  
  #note: there is no recovered node at initial state only susceptible nodes and some infected nodes
  #initial state vector must not contain 2, but only 0 and 1 for the first row
   #set next state vector initially as initial state.
  for (tick in 1:nticks) {
    inf_per_time[1,]=current_stat#time 0=initial state
    # first transition phase: S -> I(
    for (i in which(stat==0)) {#for all susceptible nodes (denoted 0) in initial state
      for (j in (1:n)) { # for all other neighbor nodes
        if ((E[i,j]>=1)&(stat[j]==1)&(runif(1)<=beta)) #each infected node infects its nearest neighbor/node it shares a link/edge with, with prob beta
        {
          nextStat[i] <- 1;break#assign node as infected if above condition is met else and break out of loop
        }
        else{
          nextStat[i]<-0 #node remains as susceptible, since not all contact leads to an infection.
        }
      }
    }
    
    # second transition phase: I -> R
    for (i in which(stat==1)) {#for all infected nodes in initial state
      if (runif(1)<=gamma) { #compares a randomly generated uniform number to recovery rate
        nextStat[i] <- 2 #assigns node as recovered if condition above is met and break out of loop else
      }
      else{
        nextStat[i]<-1 #nodes remain infected 
      }
    }

#for SIR, recovered node do not participate in further disease propagation (permanent immunity)    
    
    stat <- nextStat#assigns states of nodes for each time step as existing state 
    
    inf_per_time[time,] <- stat
    rownames(inf_per_time)=0:nticks#re-order row names from time step 0=initial state
    colnames(inf_per_time)=NULL#make column names as null
   
    
     time=time+1
    
  }
  
  return (inf_per_time)
}
G <- erdos.renyi.game(10, .5,direct=F)    #random graph generated
episim(G,beta = 0.3)


