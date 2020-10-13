library(igraph)

# (1) Function to generate a single Random graph
erdos<-function(n,p){
  erdos_graph=sample_gnp(n, p, directed = FALSE, loops = FALSE)
  summary(erdos_graph)
}


##---Test case----
f1=erdos(10,0.6)
#plot(f1)

# (2) Function for Replicating the random network N times
set.seed(1)
N=10
erdos_graph=list()
erdos<-function(n,p){
  for (i in 1:N){
    erdos_graph[[i]]=sample_gnp(n,p, directed = FALSE, loops =FALSE)
  }
  return(erdos_graph)
}