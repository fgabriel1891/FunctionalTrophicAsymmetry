
library(randomForest)
library(pROC)
###########
# Create a model to impute the multivariate trait values of non-interacting species found in JBI review

### palm traits 

# get the species score axis from the RLQ
full <- cbind(palmTraits2,RLQ$mR[,1:2])
dim(full)

# mammal traits  

# get the species score axis from the RLQ
full_m <- cbind(mammTrait5,RLQ$mQ[,1:2])
dim(full_m)


#############
## Modelling non-linear relationships with random forest
############

#### Palms 
# Impute values 
palmMiss <- palmTraits1[!rownames(palmTraits1) %in% rownames(palmTraits2),]
palmMiss2 <- palmMiss[complete.cases(palmMiss),]

dim(palmMiss2)-dim(palmMiss) # 650 species with no f.lenght or st.height info


## First axis


modelRF1 <- randomForest::randomForest(NorS1~  Stem_Height + Fruit_Length,
                                       data = full,
                                       importance =T,proximity = T, ntree = 500)


modelRF2 <- randomForest(NorS2~Stem_Height + Fruit_Length,
                         data = full)

sqrt(sum((modelRF2$predicted - full$NorS2)^2) / nrow(full))

# Number of trees that produced the lowest MSE
c(which.min(modelRF1$mse),which.min(modelRF2$mse))
sd(modelRF1$mse)
sd(modelRF2$mse)


# RMSE of the best model
c(sqrt(modelRF1$mse[which.min(modelRF1$mse)]),
  sqrt(modelRF2$mse[which.min(modelRF2$mse)]))


#### Mammals 


table(is.na(mammTrait4$Diet.Fruit))
mammMiss <- mammTrait4[!rownames(mammTrait4) %in% rownames(mammTrait5),]
mammMiss2 <- mammMiss[complete.cases(mammMiss),]
dim(mammMiss)-dim(mammMiss2)

names(full_m)

## First axis
modelRF1_m <- randomForest(NorS1~Nocturnal + Crepuscular + Diurnal + Body_Mass 
                           + Diet1 + Diet2,
                           data = full_m)
plot(modelRF1_m)
modelRF2_m <- randomForest(NorS2~Nocturnal + Crepuscular + Diurnal + Body_Mass 
                           + Diet1 + Diet2,
                           data = full_m)



# Number of trees that produced the lowest MSE
c(which.min(modelRF1_m$mse),which.min(modelRF2_m$mse))

# RMSE of the best model
c(sqrt(modelRF1_m$mse[which.min(modelRF1_m$mse)]),
  sqrt(modelRF2_m$mse[which.min(modelRF2_m$mse)]))





###############
# Make predictions 

pRF1 <- predict(modelRF1,palmMiss2)
pRF2 <- predict(modelRF2,palmMiss2)


mRF1 <- predict(modelRF1_m,mammMiss2)
mRF2 <- predict(modelRF2_m,mammMiss2)

## Extract scores and remove na

cord <- data.frame("RQL1" = mRF1, "RQL2" = mRF2)
cord$sp <- rownames(mammMiss2)
cord <- cord[complete.cases(cord),]


cord2 <- data.frame("RQL1" = pRF1,"RQL2" = pRF2)
cord2$sp <- rownames(palmMiss2 )


###
# create a new data.frame with the synthetic trait

palmStTrait <- rbind(cord2,
                     data.frame(
                       "RQL1" = full$NorS1, 
                       "RQL2"= full$NorS2, 
                       "sp" = rownames(full)))


mammalStTrait <- rbind(cord,data.frame(
  "RQL1" = full_m$NorS1, 
  "RQL2"= full_m$NorS2, 
  "sp" = rownames(full_m)))










