### Making epidemic summary on a network

# make state counts
# data frame for counts of state
# set number of rows and columns
# column names to SEIR
epi_summary<-function(f){   
  #  if (is.data.frame(f)){
  #    f=as.matrix(f)
  #  }
  numrows=dim(f)[1]
  state_counts=matrix(0,nrow=numrows,ncol=3)
  colnames(state_counts)=c("S","I","R")
  
  
  for (tick in 1:numrows){
    s = length(which(f[tick,]==0))
    i = length(which(f[tick,]==1))
    r = length(which(f[tick,]==2))
    
    state_counts[tick,] <- c(s,i,r)
  
    
  }
  
  #state_counts <- do.call(rbind, state_counts)
  state_counts <- as.data.frame(state_counts)
  rownames(state_counts)=0:(numrows-1)
  state_counts
}  



