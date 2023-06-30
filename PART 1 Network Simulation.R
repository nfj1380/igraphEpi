
#-------------------------------------------
# Calculating network properties and proportion infected from multiple graph simulations
#-------------------------------------------
#code by Nick Fountain-Jones May/June 2021

#note that to keep this succint, we only calculate for one gamma value at a time. Needs to be changed for further tests
library(igraph)
library(EpiModel)
library(intergraph)
library(DiagrammeR)
library(tidyverse)
library(doParallel)
library(foreach)
library(assortnet)

#souce functions
source('simPathE.R')
source('countStates.R ')
source('Global_net_summary.R ')
source('Network_sim_pipeline .R ')
source("estTheoreticalPrev.R")

#make multicore
cl <- makePSOCKcluster(6)
registerDoParallel(cl)

#########################################################################################
#PIPELINE
#########################################################################################

#create graph object
folder_name = c('Networks') #.edges file need to be in the working directory too at this stage.

# <- Network_sim_pipeline (folder_name)

#estimate spread across networks
summary_data2 <- Network_sim_pipeline(folder_name)

write.table(summary_data2, 'gamma0.2_summary.csv')

#print individual graphs
networks <- list.files(paste(folder_name))
filename=paste(folder_name, networks[40], sep="/") 
filename #check network name

#just get structural elements

featureSet<- getGraphFeatures(folder_name) 
write.csv(featureSet, 'FeatureSet_Jun2022')

GGally::ggpairs(featureSet)
pca_data <- na.omit(featureSet)

pca_res <- prcomp(pca_data, scale. = TRUE)

library(ggfortify)
 p <- ggplot2::autoplot(pca_res, data = pca_data, 
                        loadings = TRUE, loadings.colour = 'black',
                        loadings.label = TRUE, loadings.label.size = 3)+
  scale_colour_manual(values=c("azure3",'darkblue','green', 'deepskyblue','darkgray' ))

dat <- read.table(filename) #haven't dealt with the temporal nature of this yet. Third column appears to be a date...
g <- igraph::graph_from_data_frame(dat) %>%  as.undirected( "collapse")
l <- layout_with_fr(g)
l <- layout_with_mds(g) #a few output options

gorder(g)
library(graphlayouts)
xy <- layout_with_eigen(g,type = "adjacency",ev = "largest")

xy <- layout_with_eigen(g,type = "laplacian",ev = "smallest")
plot(g, vertex.label=NA, edge.color="orange", vertex.color="gray50", layout=l, vertex.size=8)

netclust <- clusters(g) #look for subgraphs
gcc <- V(g)[netclust$membership == which.max(netclust$csize)]#select vertices from the largest sub-graph
#we could do this for other component graphs?
gccnet <- induced.subgraph(g, gcc) #make it a igraph object.

#need to remove loops etc
gcc_s <- igraph::simplify(gccnet, remove.multiple = TRUE, remove.loops = TRUE)

plot(gcc_s , vertex.label=NA, edge.color="orange", vertex.color="gray50", layout=l, vertex.size=8)


#lots of NAs when beta is 0.01
str(test_netUpdates)
saveRDS(test_netUpdates,'test_netUpdates_gammaChanged' )

#test_gammflipped <- readRDS('test_gammflipped')
#dfNumeric <- test %>% select(-Network)

write.csv(x=summary_data , file='updated_net_char33net.csv')
#DataExplorer::create_report(dfNumeric )

##################################################################################################################
#Compare to theoretical estimates


beta = 0.1
gam = 0.04
R0 = beta/gam

Theoretical_estimates <- estTheoreticalPrev(folder_name, R0=2.5)
write.table(Theoretical_estimates, 'TheoryR0_25b.csv')

##################################################################################################################
