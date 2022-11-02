# Quantifying Functional trophic Asymmetry
# G. Muñoz
# Nov. 2022 


# load neotropics extent
neotropic <- rgdal::readOGR("DATA/Lowenberg_Neto_2014.shp", 
                            stringsAsFactors=FALSE, verbose=FALSE)


## Quantifying assymetry of interaction niches 
plot(round(PalmPresenceDatPrun$x),
     round(PalmPresenceDatPrun$y), pch = ".")


length(unique(PalmPresenceDatPrun$spName))
length(unique(MammPresenceDatPrun$spName))


# Extract trait data for mammals
mamPA_trait <- data.frame(MammPresenceDatPrun,
                          mammalStTrait[match(MammPresenceDatPrun$spName, 
                                              mammalStTrait$sp),][,1:2])

# Extract trait data for palms
palmPA_trait <- data.frame(PalmPresenceDatPrun,
                           palmStTrait[match(PalmPresenceDatPrun$spName, 
                                             palmStTrait$sp),][,1:2])

# Overlay data with neotropical provinces

mamPA_trait$Province <- over(mamPA_trait,neotropic)$Province_1
palmPA_trait$Province <- over(palmPA_trait,neotropic)$Province_1

mamPA_trait$Dominions <- over(mamPA_trait,neotropic)$Dominions
palmPA_trait$Dominions <- over(palmPA_trait,neotropic)$Dominions

# remove data outside the neotropical limits
mamPA_trait <- mamPA_trait[-which(is.na(mamPA_trait$Province)),]
palmPA_trait <- palmPA_trait[-which(is.na(palmPA_trait$Province)),]

# remove species with no trait data

palmPA_trait <- palmPA_trait[-which(is.na(palmPA_trait$RQL1 )),]
mamPA_trait <- mamPA_trait[-which(is.na(mamPA_trait$RQL1 )),]
# Add family info
mamPA_trait$Family <- mammElt$MSWFamilyLatin[match(mamPA_trait$spName,mammElt$Scientific)]
palmPA_trait$Family <- palmManNet$PalmTribe[match(palmPA_trait$spName,palmManNet$PALM)]
# remove bats as dispersers of palms 
mamPA_trait2 <- mamPA_trait[!mamPA_trait$Family == "Phyllostomidae",]
# compute richness per grid

MRich <- aggregate(mamPA_trait2$spName,list(mamPA_trait2$unID), length)
PRich <- aggregate(palmPA_trait$spName,list(palmPA_trait$unID), length)


spR <- aggregate(palmPA_trait$spName, list(palmPA_trait$unID), length)
spM <- aggregate(mamPA_trait2$spName, list(mamPA_trait2$unID), length)

###### FTA 


# Compute functional richness and functional evenness

# Functional richness
INBmgr <- sapply(unique(palmPA_trait$unID), 
                 function(x) NicheBreadth(mamPA_trait2[mamPA_trait2$unID == x,]))

INBpgr <- sapply(unique(palmPA_trait$unID), 
                 function(x) 
                   NicheBreadth(palmPA_trait[palmPA_trait$unID == x,]))

# divide it to the full extent of the pool so make it into a proportional measure

INBmgr <- INBmgr/NicheBreadth(mamPA_trait2)
INBpgr <- INBpgr/NicheBreadth(palmPA_trait)

# Functional evenness

INDmgr <- sapply(unique(palmPA_trait$unID), 
                 function(x) 
                   NicheDispersion(mamPA_trait2[mamPA_trait2$unID == x,]))

INDpgr <- sapply(unique(palmPA_trait$unID), 
                 function(x) 
                   NicheDispersion(palmPA_trait[palmPA_trait$unID == x,]))

### Visualize results

plot(INBmgr,
     INBpgr, frame = F)
abline(a=0,b=1)

plot(INDmgr,
     INDpgr, frame = F)
abline(a=0,b=1)

### Compute Z-scores to account sp. richnness in our measures of FRA


# Open parallel processing 

library(snowfall)

sfInit(parallel=TRUE, cpus=25, type="SOCK")  # start cluster 
sfSource("Scripts/FunctEcology/0_Functions.R") # Source functions to the cluster
# Export  object(s) to the cluster 
sfExport("palmPA_trait")
sfExport("spR")
sfExport("INBpgr")
sfExport("INDmgr")
sfExport("INDpgr")

sfExport("mamPA_trait2")
sfExport("spM")
sfExport("INBmgr")

ZscoresINB <- sfSapply(1:nrow(spR), function(x) makeZscore(999,palmPA_trait,spR,INBpgr,x))
ZscoresIND <- sfSapply(1:nrow(spR), function(x) makeZscore2(999,palmPA_trait,spR,INDpgr,x))

ZscoresINBrawMean <- sfSapply(1:nrow(spR), function(x) makeZscorerqw(999,palmPA_trait,spR,INBpgr,x))
ZscoresINDrawMean <- sfSapply(1:nrow(spR), function(x) makeZscoreraw2(999,palmPA_trait,spR,INDpgr,x))


ZscoresINBm <- sfSapply(1:nrow(spM), function(x) makeZscore(999,mamPA_trait2,spM,INBmgr,x))
ZscoresINDm <- sfSapply(1:nrow(spM), function(x) makeZscore2(999,mamPA_trait2,spM,INDmgr,x))

ZscoresINBmrawMean <- sfSapply(1:nrow(spM), function(x) makeZscorerqw(999,mamPA_trait2,spM,INBmgr,x))
ZscoresINDmrawMean <- sfSapply(1:nrow(spM), function(x) makeZscoreraw2(999,mamPA_trait2,spM,INDmgr,x))

sfStop()
#### close parallel processing

# Visualize results 


par(mfrow = c(1,2))
boxplot(INBpgr,ZscoresINBrawMean,INBmgr,ZscoresINBmrawMean, 
        frame = F,  lty = c(1,1,1,1), 
        col = scales::alpha(c("darkgreen", "white", "orange", "white"),0.6),
        border = c("darkgreen", "darkgreen", "orange", "orange"))
boxplot(INDpgr,ZscoresINDrawMean,INDmgr,ZscoresINDmrawMean, 
        frame = F,  lty = c(1,5,1,5), 
        ylim = c(0,1),
        col = scales::alpha(c("darkgreen", "white", "orange", "white"),0.6),
        border = c("darkgreen", "darkgreen", "orange", "orange"))

boxplot(ZscoresINB,ZscoresINBm,
        ZscoresIND,ZscoresINDm, 
        frame = F, notch = F,ylim = c(-5,10),
        col = scales::alpha(c("darkgreen","orange", "darkgreen", "orange"),0.6),
        border = c("darkgreen","orange", "darkgreen", "orange"))

dev.off()

boxplot(ZscoresINB,ZscoresINBm,
        ZscoresIND,ZscoresINDm, 
        frame = F, notch = F,ylim = c(-5,10),
        col = scales::alpha(c("darkgreen","orange", "darkgreen", "orange"),0.6),
        border = c("darkgreen","orange", "darkgreen", "orange"))

abline(h = 0)

dev.off()









