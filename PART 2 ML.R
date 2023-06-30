#install.packages("devtools")
#devtools::install_github('nfj1380/mrIML')
library(mrIML)
library(tidyverse)
library(tidymodels)
library(vip)
library(flashlight)
library(doParallel)
library(janitor)
library(viridis)
library(kernlab)
library(GGally)
library(ggfortify)
library(future.apply)
library(hrbrthemes)
library(iml)
library(plyr)

Simdata <- read.csv('NetworkDataSim_CombinedUpdated_sep22.csv')
Simdata[is.na(Simdata)] <- 0

#note that - were changed to _ to avoid downstream issues in names. 
colnames(Simdata)[2] <- c('Graph.Name')
str(Simdata)  
Graph_metaData <- read.csv('Network_repository_metaData.csv')
str(Graph_metaData)

#convert characters to factors
Graph_metaData[sapply(Graph_metaData, is.character)] <- lapply(Graph_metaData[sapply(Graph_metaData, is.character)], 
                                       as.factor)
#merge

all_data <- left_join(Simdata, Graph_metaData , by = c("Graph.Name"))

str(all_data)

#remove features that are of limited value

All_data_reducedSet <- subset(all_data, select=-c(Mammal.,Citation)) #mammal not needed as we are looking across all groups

 
#remove missing data. SOme missing responses (removes 8 networks)
All_data_reducedSetNoNA<- All_data_reducedSet[complete.cases(All_data_reducedSet  ), ] 

#remove all networks with < 10 nodes
dataNoUnder10<- All_data_reducedSetNoNA  %>% filter( Network_size >= 10) #removes 80 graphs

mean(dataNoUnder10$Network_size)

#most_infected_node removed as not a meaningful predictor in these models

Xworking <- dataNoUnder10[-c(1:21)]

Ys <- dataNoUnder10[c(2:21)]

Ys_prop <- Ys %>% select(contains("prop"))

Ys_time <- Ys %>% select(contains("time"))

Xfiltered <- subset(Xworking, select=-c(Source,  Edge_weight)) #these features aren't useful here

glimpse(Xfiltered)

dataNoCOr <- Xfiltered %>%  select(-c('Transitivity','Highest_degree', 'Adj_val','Diameter', 'Mean_degree' ))

##############################################################################################
#PCA 
##############################################################################################

pca_data <- Xworking[1:11]
pca_res <- prcomp(pca_data, scale. = TRUE)

autoplot(pca_res)

#with labels to work out outliers

autoplot(pca_res, data = Xworking, colour = 'Class',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3,
         label = TRUE)

#not labels

 p <- autoplot(pca_res, data = Xworking, colour = 'Class',
         loadings = TRUE, loadings.colour = 'black',
         loadings.label = TRUE, loadings.label.size = 3)+
          scale_colour_manual(values=c("azure3",'darkblue','green', 'deepskyblue','darkgray' ))

 p+theme_bw()

  # scale_fill_manual(colo)   

#Outlier (459) = complete voles network
#Outlier (500) ~complete tortoise network.

#kmeans

autoplot(kmeans(pca_data, 3), data = pca_data, label = TRUE, label.size = 3)

##############################################################################################
#Correlations
#############################################################################################

#species has too many levels so remove for the time being
Xfiltered <- subset(Xworking, select=-c(Source, Species, Edge_weight, Class)) #these features aren't useful here
ggpairs(Xfiltered) #lots of strong correlations


dataNoCOr <- Xfiltered %>%  select(-c('Transitivity','Diameter','Modularity', 'Mean_pathLength', 'Mean_degree', 'Class', 'Highest_degree' ))
ggpairs(dataNoCOr)
#modularity and mean degree strongly correlated for example

glimpse(dataNoCOr)
create_report(Xs)
#beta 0.5 and beta 0.1 estimates strongly correlated. Choose 0.05 


##############################################################################################
#test on a subset with the most common assocation type
##############################################################################################

#data_interactionType <- dataNoUnder10 %>% filter(Interaction_type=='Social_projection_bipartite')

data_interactionType <- dataNoUnder10 %>% filter(Interaction_type=='Physical_contact')

Xworking <- data_interactionType[-c(1:21)]

Ys <- data_interactionType[c(2:21)]

Ys_prop <- Ys %>% select(contains("prop"))

Ys_time <- Ys %>% select(contains("time"))

Xfiltered <- subset(Xworking, select=-c(Source,  Edge_weight)) #these features aren't useful here

glimpse(Xfiltered)

dataNoCOr <- Xfiltered %>%  select(-c('Transitivity','Highest_degree','Diameter', 'Mean_degree', 'Interaction_type',
                                      'Modularity', 'Class', 'Mean_pathLength'))



##############################################################################################
#ML models
##############################################################################################

model1 <- 
  rand_forest(trees = 1000, mode = "regression", mtry = tune(), min_n = tune()) %>% #100 trees are set for brevity
  set_engine("ranger", importance = "impurity") #random forest doesn't work well with dummy variables

model2 <- linear_reg() %>% 
  set_engine("lm") %>% 
  set_mode("regression")

library(kernlab)
model3<-
  svm_rbf(cost = tune(), rbf_sigma = tune()) %>%
  set_engine("kernlab") %>%
  set_mode("regression")



#X <- dataNoUnder10 %>% select(avg_prop_infected_beta0.05)

# Define set the outcomes of interest

Xraw <- dataNoCOr #up to here  - need to add species again (no -s)


#If needed - need to upadte MrIML
library(fastDummies)

Xdummy <- dummy_cols(X, remove_first_dummy = TRUE, remove_selected_columns = TRUE)

#X %>% 
 # rename(`Captive_Semi-ranging` = Captive_Semi_ranging) #this isnt working hyphen

#colnames(Ydummy[6]) <- 'Captive_Semi_ranging' #r doesn't like '-' sometimes

#set up multicore
cl <- parallel::makeCluster(5)
plan(cluster, workers=cl)

--saveRDS(yhats_rf_prop , 'yhT_prop_infected_socialBipartit_July2023a')

#save(yhats, file='logreg_model')
ModelPerf <- mrIMLperformance(yhats_rf_prop, model1, Y=Ys_prop) 
ModelPerf 

#extract predictions
predictionsR05 <- collect_predictions(yhats_rf_prop[[1]]$last_mod_fit, summarize = TRUE, , grid[1, ]) %>% arrange(.row)

training_pred <- 
  predict(yhats_rf_prop[[1]]$mod1_k, as.data.frame(yhats_rf_prop[[1]]$data_train)) %>% 
  bind_cols(predict(yhats_rf_prop[[1]]$mod1_k,yhats_rf_prop[[1]]$data_train)) %>% 
  # Add the true outcome data back in
  bind_cols(as.data.frame(yhats_rf_prop[[1]]$data_train) %>% 
              select(class))

#predicts well on test data
testing_pred <- 
  predict(yhats_rf_prop[[1]]$mod1_k, as.data.frame(yhats_rf_prop[[1]]$data_testa)) %>% 
  bind_cols(predict(yhats_rf_prop[[1]]$mod1_k,yhats_rf_prop[[1]]$data_testa)) %>% 
  # Add the true outcome data back in
  bind_cols(as.data.frame(yhats_rf_prop[[1]]$data_testa) %>% 
              select(class))

#should now be fixed - work off mrvip_v2
VI <- mrVip (yhats_rf_prop, X=Xraw)

 p <- plot_vi(VI=VI,  X=Xraw,Y=Ys_prop, modelPerf=ModelPerf, cutoff= 0)+theme_bw()
 #checks which particular species etc
 
flashlightObj_rf <- mrFlashlight(yhats_rf_prop, X=Xraw,Y=Ys_prop, response = "multi", mode='regression')

#Interpretation

aledata_spec <-light_profile(flashlightObj_rf, v = "Spectral_radius", type = "ale",  n_bins =50)

mrProfileplot(aledata_spec , sdthresh = 0)+
  theme_bw()+
  geom_rug(data=Xraw, aes(x =Spectral_radius), inherit.aes = F)


aledata_Fiedler <-light_profile(flashlightObj_rf, v = "Fiedler", type = "ale",  n_bins =50)

mrProfileplot(aledata_Fiedler , sdthresh = 0)+
  theme_bw()+
  geom_rug(data=Xraw, aes(x =Fiedler), inherit.aes = F)

aledata_Cent <-light_profile(flashlightObj_rf, v = "Centrality", type = "ale",  n_bins =50)

mrProfileplot(aledata_Cent , sdthresh = 0)+
  theme_bw()+
  geom_rug(data=X, aes(x =Centrality), inherit.aes = F)

aledata_Mod <-light_profile(flashlightObj_rf, v = "Qrel", type = "ale",  n_bins =50)

mrProfileplot(aledata_Mod , sdthresh = 0)+
  theme_bw()+
  geom_rug(data=X, aes(x =Modularity), inherit.aes = F)

#nteractions

#cl <- parallel::makeCluster(5)
#plan(cluster, workers=cl)

#seemingly only work after running MrIML predicts
interactions <-mrInteractions(yhats_rf_prop, Y=Ys_prop,X=X,  mode='regression') #this is computationally intensive so multicores are needed.


mrPlot_interactions(interactions, Y=Ys_prop,X=X, top_ranking = 5, top_response=10)

mrIMLconverts_list <- MrIMLconverts(yhats_rf_prop,X=X, mode='regression')

featureA = 'Centrality' #doesn't do catergories
featureB = 'Fiedler'

p2d <- mrProfile2D (mrIMLconverts_list, featureA,
                    featureB,  mode='regression',
                    grid.size=30, method = 'ale')
plot(p2d) +theme_bw()

str(X)

#plot prediction scatter for all responses.

plot(light_profile(flashlightObj, v = "Modularity", type = "ale",  n_bins =25))+
  theme_bw()+
  geom_rug(data=Y, aes(x =Modularity), inherit.aes = F)

plot(light_profile(flashlightObj, v = "network_size", type = "ale", n_bins =100))+
  theme_bw()+
  geom_rug(data=Y, aes(x =network_size), inherit.aes = F)

plot(light_profile(flashlightObj, v = "Centrality", type = "ale",  n_bins =50))+
  theme_bw()+
  geom_rug(data=Y, aes(x =Centrality), inherit.aes = F)
  
plot(light_profile(flashlightObj, v = "Fiedler", type = "ale", n_bins =100))+
  theme_bw()+
  geom_rug(data=Y, aes(x=Fiedler), inherit.aes = F)

plot(light_ice(flashlightObj, v = "Fiedler", center = "first"))+
  theme_bw()+
  geom_rug(data=Y, aes(x =Fiedler), inherit.aes = F)

plot(light_ice(flashlightObj, v = "Centrality", center = "first"))+
  theme_bw()+
  geom_rug(data=Y, aes(x = Centrality), inherit.aes = F)

plot(light_profile2d(flashlightObj, v = c("network_size", "Centrality")))

plot(light_profile2d(flashlightObj, v = c("network_size", "Mean_degree")))
geom_rug(data=Y, aes(x = network_size, y=Mean_degree), inherit.aes = F)+
  scale_fill_continuous(high = "#56B1F7", low = "#132B43")

plot(light_profile2d(flashlightObj, v = c("network_size", "Fiedler")))+
  geom_rug(data=Y, aes(x = network_size, y=Fiedler), inherit.aes = F)+
  #scale_fill_gradientn(colours = colorspace::heat_hcl(7))
  scale_fill_continuous(high = "#132B43", low = "#56B1F7")

plot(light_global_surrogate(flashlightObj))

st <- light_interaction(flashlightObj, grid_size = 30, n_max = 50, seed = 42, pairwise=T)+
  theme_bw()
plot(st)
