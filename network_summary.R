#[2] simulating SIR process on each simulated graph and extracting the summary

sim_all_graphs=NULL
summary_all_graphs=NULL

nets_epi_summary<-function(G){
  for (i in 1:N){
    sim_all_graphs[[i]] = lapply(G,episim)
    summary_all_graphs[[i]]=lapply(sim_all_graphs[[i]],epi_summary)
  }
  return(sim_all_graphs)
  #  return(summary_all_graphs)
}


#[3] Time to infection of 20% of individual
