
# How well they cover the metaweb with the full web

palmTraits2
mNBat <- mammTrait[!mammTrait$Family == "PHYLLOSTOMIDAE",]
mNBat <- mNBat[!mNBat$Diet.Fruit == 0,]
ProvDistMammal$Family <- mammTrait$Family[match(ProvDistMammal$L1,
                                                mammTrait$MammalName)]
MaNeo <- ProvDistMammal[!ProvDistMammal$Family == "PHYLLOSTOMIDAE",]

##### Plot trait equivalence between full and subset datasets 

par(mfrow = c(2,2), mar = c(4,2,2,2))
## Mammals 
plot(density(na.omit(log1p(mNBat$Mass))), 
     ylim = c(0,0.4), frame = F, 
     main = "", xlab = "Body mass (log)")
lines(density(
  na.omit(log(mNBat$Mass[mNBat$MammalName 
                         %in% unique(MaNeo$L1)]))),
  col = "red")
legend("topright",
       bty = "n",
       cex = 0.6,
       lty = c(1,1), col = c("red", "black"),
       c("subset wt interactions", "neotropics"))



plot(density(na.omit((mNBat$Diet.Fruit))), 
     ylim = c(0,0.1), frame = F, 
     main = "", xlab = " % frugivore")
lines(density(
  na.omit((mNBat$Diet.Fruit[mNBat$MammalName 
                            %in% unique(MaNeo$L1)]))),
  col = "red")
legend("topright",
       bty = "n",       cex = 0.6,
       lty = c(1,1), col = c("red", "black"),
       c("subset wt interactions", "neotropics"))


## Palms 

ProvDistPalm$Family <- palmTraits$PalmSubfamily[match(ProvDistPalm$L1,
                                                      palmTraits$SpecName)]

plot(density(na.omit(log1p(ProvDistPalm$FL))), 
     ylim = c(0,2), frame = F, 
     main = "", xlab = "Fruit Lenght (log)")
lines(density(
  na.omit(log1p(ProvDistPalm$FL[ProvDistPalm$L1 
                                %in% rownames(palmTraits2)]))),
  col = "red")
legend("topright",
       bty = "n",
       cex = 0.6,
       
       lty = c(1,1), col = c("red", "black"),
       c("subset wt interactions", "neotropics"))

plot(density(na.omit(log1p(ProvDistPalm$SH))), 
     ylim = c(0,1), frame = F, 
     main = "", xlab = "Stem height (log)")
lines(density(
  na.omit(log1p(ProvDistPalm$SH[ProvDistPalm$L1 
                                %in% rownames(palmTraits2)]))),
  col = "red")
legend("topright",
       bty = "n",
       cex = 0.6,
       lty = c(1,1), col = c("red", "black"),
       c("subset wt interactions", "neotropics"))
##########

dim(N)
bipartite::networklevel(N, index = "binary")
mean(rowSums(N))
mean(colSums(N))



Ng <- cassandRa::CreateListObject(nRLQ)
dsd <- cassandRa::FitAllModels(Ng)

#### 
# test accuracy with Youden's J
#### 



TestYJ <- function(probNet,obs, n){
  sq <- seq(range(probNet)[1],
            range(probNet)[2], 
            diff(range(probNet))/n)
  sens <- c()
  speci <- c()
  YJ <- c()
  
  for(i in 1:n){
    prob10 <- ifelse(probNet>sq[i], 1,0)
    
    Ttab <- prop.table(table(obs,prob10))
    
    sens[i] <- Ttab[4]/c(Ttab[4] + Ttab[2])
    speci[i] <- Ttab[1]/c(Ttab[1] + Ttab[3])
    
    YJ[i] <- sens[i]+speci[i]-1
    
  }
  
  return(data.frame(sens, speci, YJ))
  
}


## pull RQL probmat and reconstitute
head(prepNet)

hist(prepNet$value)
RLQ$mR
RLQ$mQ
prepNet <- reshape::melt(dsd$obs)
prepNet <- data.frame(prepNet,data.frame(RLQ$mQ[match(prepNet$FRUGIVORE,
                                                      rownames(RLQ$mQ)),][,1:2],
                                         RLQ$mR[match(prepNet$PALM, rownames(RLQ$mR)),][,1:2]))

prepNet$dist <- sqrt(c(prepNet$NorS1-prepNet$NorS1.1)^2 
                     + c(prepNet$NorS2-prepNet$NorS2.1)^2)
prepNet$intNet <- exp(-prepNet$dist/(4*0.07^2))

probNetRQL <- xtabs(intNet~PALM + FRUGIVORE, prepNet)



probNet <- list(dsd$SBM_ProbsMat, 
                dsd$C_ProbsMatrix, 
                dsd$M_ProbsMatrix, 
                dsd$B_ProbsMat)

YJtestin <- lapply(1:4, function(i) TestYJ(probNet[[i]],dsd$obs, 100))



plot(YJtestin[[1]]$YJ~seq(range(probNet[[1]])[1],
                          range(probNet[[1]])[2], 
                          diff(range(probNet[[1]]))/100)[-101], type = "l")




plot(YJtestin[[1]]$sens~
       c(1-YJtestin[[1]]$speci), 
     type = "l",
     frame = F,
     xlim = c(0,1),
     ylim = c(0,1))
points(YJtestin[[2]]$sens~
         c(1-YJtestin[[2]]$speci), 
       type = "l", 
       lty = 2)
points(YJtestin[[3]]$sens~
         c(1-YJtestin[[3]]$speci),
       type = "l", 
       lty = 3)
points(YJtestin[[4]]$sens~
         c(1-YJtestin[[4]]$speci),
       type = "l", lty = 4)
abline(a=0,b=1, lty = 5, col = "gray")
legend("bottomright", lty = 1:5,c("SBM", "C", "M", "MC", "Random"), 
       col = c(rep("black",4), "gray"))


####

plot(YJtestin[[1]]$YJ~seq(range(probNet[[1]])[1],
                          range(probNet[[1]])[2], 
                          diff(range(probNet[[1]]))/100)[-101], 
     type = "l", 
     ylab= "Youden's J",
     xlab = "threshold")
points(YJtestin[[2]]$YJ~seq(range(probNet[[2]])[1],
                            range(probNet[[2]])[2], 
                            diff(range(probNet[[2]]))/100)[-101], 
       type = "l", lty=2)

points(YJtestin[[3]]$YJ~seq(range(probNet[[3]])[1],
                            range(probNet[[3]])[2], 
                            diff(range(probNet[[3]]))/100)[-101], 
       type = "l",lty=3)
points(YJtestin[[4]]$YJ~seq(range(probNet[[4]])[1],
                            range(probNet[[4]])[2], 
                            diff(range(probNet[[4]]))/100)[-101], 
       type = "l", lty=4)
legend("topright", lty = 1:4,c("SBM", "C", "M", "MC"), 
       col = c(rep("black",4)))

##################
# Fit latent models 

# Get SBM

SBMs <- cassandRa::FitSBM(Ng)

## 

PalmNet <- data.frame(Ng$HostNames,
                      "SBMs.SB_H" =SBMs$SBM1$SB_H)

MammNet <- data.frame(Ng$WaspNames,
                      "SBMs.SB_W" = SBMs$SBM1$SB_W)


### Inferring the position with traits 

PalmNet <- data.frame(PalmNet,
                      palmTraits[match(PalmNet$Ng.HostNames,
                                       palmTraits$SpecName),])
MammNet <- data.frame(MammNet,
                      mammTrait[match(MammNet$Ng.WaspNames,
                                      mammTrait$Scientific),])

heatmap(table(PalmNet$accGenus,
              (PalmNet$SBMs.SB_H)), 
        col = viridis::cividis(100))

heatmap(table(MammNet$Family,
              (MammNet$SBMs.SB_W)),
        col = viridis::cividis(100))



### Predict 
names(MammNet)
MammNet$ForStrat.Value

heatmap(log1p(table(PalmNet$SBMs.SB_H,PalmNet$PalmTribe)))


medBM <- aggregate(MammNet$BodyMass.Value, list(MammNet$Family), median)
ordm <- medBM$Group.1[order(medBM$x)]

medPS <- aggregate(PalmNet$MaxFruitLength_cm, list(PalmNet$PalmTribe), median)
ordmPS <- medPS$Group.1[order(medPS$x)]

mmAso <- table(MammNet$Family,MammNet$SBMs.SB_W)
PmmAso <- table(PalmNet$PalmTribe,PalmNet$SBMs.SB_H)

heatmap(prop.table(mmAso)[match(ordm,rownames(mmAso)),],Rowv = NA,Colv = NA)
heatmap(prop.table(PmmAso)[match(ordmPS,rownames(PmmAso)),],Rowv = NA,Colv = NA)




mod1 <- nnet::multinom(SBMs.SB_W~Family,
                       data = MammNet)
mod2 <- nnet::multinom(SBMs.SB_W~Family+
                         log1p(BodyMass.Value),
                       data = MammNet)
mod3 <- nnet::multinom(SBMs.SB_W~log1p(BodyMass.Value),
                       data = MammNet)
mod4 <- nnet::multinom(SBMs.SB_W~Diet.Fruit,
                       data = MammNet)

mod5 <- nnet::multinom(SBMs.SB_W~Family+Diet.Fruit+ForStrat.Value+
                         log1p(BodyMass.Value),
                       data = MammNet)



# evaluate models 

AICdist <- (c(mod1$AIC,mod2$AIC,mod3$AIC, mod4$AIC, mod5$AIC))

modlist <- list(mod1, mod2, mod3, mod4, mod5)
# 1-diagonal of the prop.table 
EvalMod <- data.frame(AICdist,
                      "Diag" = sapply(1:5, 
                                      function(x)
                                        c(1- sum(diag(
                                          prop.table(
                                            table(MammNet$SBMs.SB_W,
                                                  predict(modlist[[x]], 
                                                          MammNet)))))))
)

plot(EvalMod)

EvalMod[order(EvalMod$Diag, decreasing = F),]

heatmap(table(MammNet$SBMs.SB_W,
              predict(modlist[[5]], 
                      MammNet)), Rowv = NA, Colv = NA)

#### Predicting palms 

names(PalmNet)
Pmod1 <- nnet::multinom(SBMs.SB_H~accGenus,data = PalmNet)
Pmod2 <- nnet::multinom(SBMs.SB_H~accGenus+
                          (AverageFruitLength_cm),
                        data = PalmNet)
Pmod3 <- nnet::multinom(SBMs.SB_H~(AverageFruitLength_cm),
                        data = PalmNet)
Pmod4 <- nnet::multinom(SBMs.SB_H~(AverageFruitLength_cm)+
                          MaxStemHeight_m,
                        data = PalmNet)

Pmod5 <- nnet::multinom(SBMs.SB_H~accGenus+MaxStemHeight_m+
                          (AverageFruitLength_cm) + Acaulescent+Erect+
                          Climbing,
                        data = PalmNet)

Pmod6 <- nnet::multinom(SBMs.SB_H~accGenus+MaxStemHeight_m+
                          (AverageFruitLength_cm) + UnderstoreyCanopy,
                        data = PalmNet)

names(PalmNet)

SBMs$SBM1$Omega_rs
# evaluate models 

AICdist_p <- (c(Pmod1$AIC,Pmod2$AIC,Pmod3$AIC, Pmod4$AIC, Pmod5$AIC,
                Pmod5$AIC))

modlist_p <- list(Pmod1, Pmod2, Pmod3, Pmod4, Pmod5, Pmod6)
# 1-diagonal of the prop.table 
EvalMod_p <- data.frame(AICdist_p,
                        "Diag" = sapply(1:6, 
                                        function(x)
                                          1- sum(diag(
                                            prop.table(
                                              table(PalmNet$SBMs.SB_H,
                                                    predict(modlist_p[[x]], 
                                                            PalmNet))))))
)


1- sum(diag(
  prop.table(
    table(PalmNet$SBMs.SB_H,
          predict(modlist_p[[1]], 
                  PalmNet))
  )
))


plot(EvalMod_p)
EvalMod_p[order(EvalMod_p$Diag, decreasing = F),]

heatmap(table(PalmNet$SBMs.SB_H,
              predict(modlist_p[[6]], 
                      PalmNet)), Rowv = NA, Colv = NA)




coef(Pmod6)[,-1]

pred.probs <- predict(Pmod6, type = "probs")
pred.probs2 <- predict(mod5, type = "probs")


plot(pred.probs[1,], type = 'l')




dt <-data.frame("Tribe" = as.numeric(as.factor(PalmNet$PalmTribe)), 
                "Stem" = PalmNet$MaxStemHeight_m, 
                "Fruit" = PalmNet$AverageFruitLength_cm,
                "Strata" = as.numeric(as.factor(PalmNet$UnderstoreyCanopy)))

rownames(dt) <- PalmNet$SpecName
dt2 <-data.frame("Family" = as.numeric(as.factor(MammNet$Family)), 
                 "Frugivory %" = MammNet$Diet.Fruit, 
                 "Body Size" = log(MammNet$BodyMass.Value),
                 "Strata" = as.numeric(as.factor(MammNet$ForStrat.Value)))


plot(vegan::rda( pred.probs), 
     cex = .4, col = "gray",
     type = "t")


plot(vegan::rda(dt2, pred.probs2), 
     cex = .4, col = "gray",
     type = "t")
str(data.frame(as.factor(PalmNet$PalmTribe), 
               PalmNet$MaxStemHeight_m, 
               PalmNet$AverageFruitLength_cm,
               as.factor(PalmNet$UnderstoreyCanopy)))
##### 
# expand 
palmTraits_test <-palmTraits[palmTraits$PalmTribe %in% PalmNet$PalmTribe,]
palmTraits_test <-palmTraits_test[palmTraits_test$UnderstoreyCanopy
                                  %in% PalmNet$UnderstoreyCanopy,]

PalmPreds <- data.frame("spNamePalm" = palmTraits_test$SpecName,
                        "group" = predict(modlist_p[[6]],
                                          palmTraits_test))



# predict mammal
mammTrait[mammTrait$Family %in% MammNet$Family,]
mammTrait_test <- mammTrait_test[mammTrait_test$ForStrat.Value
                                 %in% MammNet$ForStrat.Value,]
#
mammPreds <- data.frame("spNameMam" = mammTrait_test$Scientific,
                        "group" = predict(modlist[[5]],
                                          mammTrait_test))


####
# project geographically palms
PalmPresenceDat$group <- PalmPreds$group[match(PalmPresenceDat$spName,
                                               PalmPreds$spNamePalm)]
PalmPresenceDat_test<-PalmPresenceDat[!is.na(PalmPresenceDat$group),]
PalmPresenceDat_test$ID2 <- paste(PalmPresenceDat_test$unID,
                                  PalmPresenceDat_test$spName, sep = "_")
# summarize dataset palms
PalmPresenceDat_test <- PalmPresenceDat_test[!duplicated(PalmPresenceDat_test$ID2),]

MatrixPalm <- xtabs(~group +unID,PalmPresenceDat_test)
MatrixPalm <- t(MatrixPalm)
MatrixPalm <- MatrixPalm

matDist <- dist(MatrixPalm)
dim(matDist)

####
# project geographically mammals
MammPresenceDat$group <- mammPreds$group[match(MammPresenceDat$spName,
                                               mammPreds$spNameMam)]

MammPresenceDat_test<-MammPresenceDat[!is.na(MammPresenceDat$group),]

MammPresenceDat_test$ID2 <- paste(MammPresenceDat_test$unID,
                                  MammPresenceDat_test$spName, sep = "_")
# summarize dataset palms
MammPresenceDat_test <- MammPresenceDat_test[!duplicated(MammPresenceDat_test$ID2),]

MatrixMamm <- xtabs(~group +unID,MammPresenceDat_test)
MatrixMamm <- t(MatrixMamm)
MatrixMamm <- MatrixMamm
MatrixMammP <- MatrixMamm[rownames(MatrixMamm) %in% 
                            rownames(MatrixPalm),]


##### Reconstitute probabilistic networks

MammPresenceDat$type <- rep("M", length( MammPresenceDat$x))
PalmPresenceDat$type <- rep("P", length( PalmPresenceDat$x))

PresAll <- rbind(MammPresenceDat,PalmPresenceDat)
head(PresAll)
PresAll$ID2 <- paste(PresAll$unID, PresAll$spName, sep = "_")
# remove duplicates
PresAll <- PresAll[!duplicated(PresAll$ID2),]
dim(PresAll)
# remove NA groups 
PresAll <- PresAll[-which(is.na(PresAll$group)),]
head(PresAll)

### recover metaweb

N1 <- PresAll_prunned
M <- N1$spName[N1$type == "M"]
P <- N1$spName[N1$type == "P"]
M_P <- expand.grid(unique(M),unique(P))
M_P$G1 <- N1$group[match(M_P$Var1, N1$spName)]
M_P$G2 <- N1$group[match(M_P$Var2, N1$spName)]
# compute int probabilities 
M_P$intPro <- sapply(1:length(M_P$Var1), function(i) 
  (SBMs$SBM1$Omega_rs[M_P$G1[i], M_P$G2[i]]))
M_P$intPro2 <- abs(rnorm(length(M_P$intPro), M_P$intPro, 0.03))
FullnetT <- xtabs(intPro2~Var1 + Var2, M_P)




#######
points(-41,-10, col = "red", lwd = 3) #FEA_t
points(-59,6, col = "red", lwd = 3) # FEA_f
points(-53,5, col = "red", lwd = 3) # FRA_f
points(-57,-25, col = "red", lwd = 3) # FRA_t

head(rownames(trySEM)[order(trySEM$H2, decreasing = T)])

plot(trySEM$H2[order(trySEM$H2, decreasing = F)])
Nnetsub <- M_P[M_P$Var2 %in% unique(PalmPresenceDatPrun$spName[PalmPresenceDatPrun$unID == "-72_-12"])&
                 M_P$Var1 %in% unique(MammPresenceDatPrun$spName[MammPresenceDatPrun$unID == "-72_-12"]),]


FullnetT2 <- xtabs(intPro~Var1 + Var2, droplevels(Nnetsub))

head(FullnetT)
dim(FullnetT2)


(fullNetGraph)
fullNetGraph <- igraph::graph_from_incidence_matrix(FullnetT2)

lm <- length(c(na.omit(mammalStTrait$RQL1[match(names(V(fullNetGraph)),mammalStTrait$sp)])))
lp <- length(na.omit(palmStTrait$RQL1[match(names(V(fullNetGraph)),palmStTrait$sp)]))


coordsNet <- data.frame("R1" = c(na.omit(mammalStTrait$RQL1[match(names(V(fullNetGraph)),mammalStTrait$sp)]),
                                 na.omit(palmStTrait$RQL1[match(names(V(fullNetGraph)),palmStTrait$sp)])),
                        "R2"= c(na.omit(mammalStTrait$RQL2[match(names(V(fullNetGraph)),mammalStTrait$sp)]),
                                na.omit(palmStTrait$RQL2[match(names(V(fullNetGraph)),palmStTrait$sp)])))


ew <- M_P$intPro[match(attr(E(fullNetGraph), "vnames"),
                       paste0(M_P$Var1,"|",M_P$Var2))]




plot(fullNetGraph, 
     vertex.size =5, 
     vertex.label = "",
     vertex.color = c(rep("red",lm),
                      rep("green", lp)),
     edge.color = f(ew, 4,"Spectral"),
     edge.width = (ew),
     layout = as.matrix(coordsNet))
hist(ew)







#####
# Recreate functional domain
######

### P

########
# recover networks per grid
#########

rowSums(as.matrix(SBMs$SBM1$Omega_rs))
apply(SBMs$SBM1$Omega_rs, 1, sd)
sd(SBMs$SBM1$Omega_rs)

sd(diag(SBMs$SBM1$Omega_rs))
mean(diag(SBMs$SBM1$Omega_rs))

intMat <- SBMs$SBM1$Omega_rs
diag(intMat) <- NA
sd(intMat, na.rm = T)
mean(table(SBMs$SBM1$SB_W))
barplot(rbind(c(table(SBMs$SBM1$SB_H)),c(table(SBMs$SBM1$SB_W))), 
        col = c("#70AD47", "#BE5107"))

plot(sort(SBMs$SBM1$SB_H), type = "l", xlim = c(0,100))
points(sort(SBMs$SBM1$SB_W), type = "l")




graphSBM <- igraph::graph_from_adjacency_matrix(as.matrix(SBMs$SBM1$Omega_rs), 
                                                mode = "undirected" ,
                                                weighted = T)

l <- layout_with_fr(graphSBM, weights=E(graphSBM)$weight)
plot(graphSBM, layout = l, edge.width = E(graphSBM)$weight*20)



CalCMetric <- function(PresAll_prunned,SBMs, sdErr,x, rarSam){
  N1 <- PresAll_prunned[PresAll_prunned$unID == 
                          unique(PresAll_prunned$unID)[x],]
  
  N2 <- expand.grid(N1$spName[N1$type == "M"],
                    N1$spName[N1$type == "P"])
  
  print(tail(N1))
  N2$G1 <- N1$group[match(N2$Var1, N1$spName)]
  N2$G2 <- N1$group[match(N2$Var2, N1$spName)]
  # compute int probabilities 
  N2$intPro <- sapply(1:length(N2$Var1), function(i) 
    (SBMs$SBM1$Omega_rs[N2$G1[i], N2$G2[i]]))
  
  N2$intPro2 <- abs(rnorm(length(N2$intPro), N2$intPro, 0.03))
  netT <- xtabs(intPro2~Var1 + Var2, N2)
  netMet <- cassandRa::RarefyNetwork(netT,
                                     abs_sample_levels = rarSam,
                                     metrics = c("SA", "connectance",
                                                 "web asymmetry", 
                                                 "ISA", "H2"),
                                     ... = list(H2_integer = F),
                                     output = "CI")
  return(netMet)
  
}

palmTab <- table(PalmPresenceDat$spName, PalmPresenceDat$unID)
palmTab[palmTab>1]<-1
mamlTab <- table(MammPresenceDatPrun$spName, MammPresenceDatPrun$unID)
mamlTab[mamlTab>1]<-1

RMam <- colSums(mamlTab)
RPal <- colSums(palmTab)

RAss <- RPal[match(names(RMam), names(RPal))]-RMam
names(RAss) <- names(RMam)

RAss <- na.omit(RAss)


### check palm and mammal p/a data 


PMtb <-as.data.frame.matrix(t(table(PresAll$type, PresAll$unID)))
# remain only with those coordinates which have at least 5 palms and 5 mammals 
PMtb <- PMtb[PMtb$P > 5 & PMtb$M > 5,]

PresAll_prunned <- PresAll[PresAll$unID %in% rownames(PMtb),]


#########
# calculate network metrics for a sample of grids
######


median(table(PresAll$unID[PresAll$type == "P"]))
median(table(PresAll$unID[PresAll$type == "M"]))


## computation needs cluster to be effective 

library(snowfall)

sfInit(parallel=TRUE, cpus=30)
sfExport("PresAll_prunned") # export object to cluster 
sfExport("SBMs") # export object to cluster 
sfSource("Scripts/FunctEcology/016_FinalManuscript/ReviewResubmit/00_functionsExtra.R")

RarefMetrics <- sfLapply(1:length(unique(PresAll_prunned$unID)), 
                         function(i) 
                           CalCMetric(PresAll_prunned, SBMs, 0.03, i, 100))

sfStop() ## stop cluster 

names(RarefMetrics) <- unique(PresAll_prunned$unID)
# turn list into dataframe
Rmetrics <- reshape2::melt(RarefMetrics)
# get only the means of each metric
Rmetrics_mean <- Rmetrics[Rmetrics$variable == "Mean",]

# add lat_long 
Rmetrics_mean <- data.frame(Rmetrics_mean,
                            apply(
                              stringr::str_split(Rmetrics_mean$L1, 
                                                 "_", simplify = T),
                              2, as.numeric))

unique(Rmetrics_mean$Metric)

SA_mean <- Rmetrics_mean[Rmetrics_mean$Metric == "specialisation asymmetry",]
H2_mean <- Rmetrics_mean[Rmetrics_mean$Metric == "H2",]
Co_mean <- Rmetrics_mean[Rmetrics_mean$Metric == "connectance",]


min(H2_mean$value)

##################################

