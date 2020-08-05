##

getGraphFeatures <- function() {
  ## output: data frame
  
}

# converting downloaded networks into igraph objects
# use something liek this to convert:
#dat <- read.table(networks[i]) #haven't dealt with the temporal nature of this yet. Third column appears to be a date...
#g <- graph_from_data_frame(dat)%>% as.undirected( "collapse")

G <- graph_from_data_frame(dat)%>% as.undirected( "collapse")


#Simple SIR disease simulation
# Disease status: 0=healthy/susceptible, 1=infected, 2=removed/dead/recovered
episim <- function(G, nticks=300, beta=0.1, gamma=0.2, initState=NULL, propInfected=0.5) {
  ## G is igraph object
  E <- G[,] #adjacency matrix of graph (network object)
  n <- ncol(E) #number of columns (nodes) in network object
  
  # Set initial state  
  # initState is list of initial state of each vertex in G.  E.g., from (1,0,0...0) or (0,1,0,1,0,....)
  if (is.null(initState)) {
    stat <- rbinom(n,1,propInfected)#set initial state of each nodes if non (note, no recovered node at initial state)
  } else {
    stat <- initState #initial existing state
  }
  
  nextStat <- rep(0,length(stat))#set vector to length of initial state to store node states in each step
  
  for (tick in 1:nticks) {
    
    # first transition phase: S -> I(
    for (i in which(stat==0)) {#for all susceptible nodes (denoted 0) in initial state
      for (j in (1:n)) { # for all other neighbor nodes
        if ((E[i,j]>=1)&(stat[j]==1)&(runif(1)>=beta)) #each infected node infects its nearest neighbor/node it shares a link/edge with, with prob beta
        {
          nextStat[i] <- 1 #assign node as infected if above condition is met else
        }
        else{
          nextStat[i]<-0;break #node remains as susceptible, since not all contact leads to an infection, and break out of loop
        }
      }
    }
    
    # second transition phase: I -> R
    for (i in which(stat==1)) {#for all infected nodes in initial state
      if (runif(1)>=gamma) {#compares a randomly generated uniform number to recovery rate
        nextStat[i] <- 2 #assigns node as recovered if condition above is met else
      else{
        nextStat[i]<-1;break #nodes remain infected and breaks out of loop
      }
      }
    }
#for SIR, recovered node do not participate in further disease propagation (permanent immunity)    
    
    stat <- nextStat#assigns states of nodes for each time step as existing state 
    #print(length(which(stat==1)))
    print(stat)
  }
}

G <- erdos.renyi.game(20, .10,direct=F)    #random graph generated

episim(G)


