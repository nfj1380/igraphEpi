library(igraph)

#(1) Function to generate a single scale-free graph
scale_free<-function(n,pw){
  scale_free_graph=sample_pa(n,pw, directed= F)
  summary(scale_free_graph)
}


#-----Test case-----
f4=scale_free(50,1)
#plot(f4)



# (2) Function for Replicating the scale-free network N times

set.seed(1)
N=10
scale_free_net=list()
scale_free<-function(n,pw){
  for (i in 1:N){
    scale_free_net[[i]]=sample_pa(n,pw, directed= F)
  }
  return(scale_free_net)
}

#--Test case----
f4=scale_free(50,1)
#plot(f4)