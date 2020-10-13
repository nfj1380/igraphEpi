#[1] simulating a single SIR process on each simulated graph and extracting the summary

single_epi_summary<-function(f){
    summary_all_graphs=lapply(f,epi_summary)
}


#[2] simulating multiple SIR processes on each simulated graph and extracting the summary

summary_all_graphs=NULL
N=10

nets_epi_summary<-function(f){
    for (i in 1:N){
      summary_all_graphs[[i]]=lapply(f,epi_summary)
    }
    return(summary_all_graphs)
}



#[3] Time to infection of 20% of individual
