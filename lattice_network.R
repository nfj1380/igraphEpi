library(igraph)

# (1) Function to generate a single lattice graph
Lattice<-function(d,l){
  lattice_graph=make_lattice(dim=d,length=l, directed= F, mutual= T, circular = F)
  summary(lattice_graph)
}



# (2) Function for Replicating the lattice network N times
N=10
lattice_graph=list()
Lattice<-function(d,l){
  for (i in 1:N){
    lattice_graph[[i]]=make_lattice(dim=d,length=l, directed= F, mutual= T, circular = F)
  }
  return(lattice_graph)
}

##---Test case----
f3=Lattice(2,5)
#plot(f3)