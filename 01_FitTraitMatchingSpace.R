## Fit a trait matching space 
# Gabriel Muñoz 
# Nov 2022
## Load data 
# Palm-frugivore interactions 
palmMetaNet <- read.csv("DATA/PalmDryadRepo/PalmFrugDatasetOCT2018.csv", 
                        header = T, stringsAsFactors = F)

# bird traits
birdsTr <- read.csv("DATA/VAR/BIRDS/BirdFuncDat.txt",
                    sep = "\t",
                    header = T,
                    stringsAsFactors = F)
# palm traits
palmTraits <- read.csv("DATA/PalmTraits_1.0.txt", sep = "\t", stringsAsFactors = F)

# mammal traits
# body mass
mammMass <- read.csv("DATA/VAR/Mammals/Mammal_BodyMass/MammalMassSandom2013.csv", 
                     sep = ",",
                     header = T,
                     stringsAsFactors = F)
# elton traits 
mammElt <- read.csv("DATA/VAR/Mammals/Elton Traits/MamFuncDat.txt",
                    sep = "\t",
                    header = T,
                    stringsAsFactors = F)


############
## Generating multidimensional trait space RQL axis (based on interactions) adapted from J.A script
############

################################################################################
# load packages and data
################################################################################
library(bipartite)
library(ade4)
library(FD)
library(picante)



# erase unwnated columns 
names(palmTraits1)
palmTraits1 <- palmTraits[-c(1,2,3,4,5,29,25,26,27,28)]

palmTraits1 <- palmTraits[c( "MaxStemHeight_m", "AverageFruitLength_cm")]
rownames(palmTraits1) <- palmTraits$SpecName

dim(palmTraits1)
head(palmTraits1)
# trasnform to cm
palmTraits1$MaxStemHeight_m <- palmTraits1$MaxStemHeight_m*100

# standardization and normal transformation 

colnames(palmTraits1) <- c("Stem_Height",  "Fruit_Length")


### Mammals 

# match names between mammal datasets 
mammElt <- mammElt[match(mammMass$MammalName,mammElt$Scientific),]
# bind mammal trait information 
mammTrait <- cbind(mammElt, mammMass)
# erase unwanted columns
mammTrait1 <- mammTrait[,-c(3,29,22, 16,15,30,14, 15,17,18,22,23,31,25,26,27,28,1,2,3)]
# give names
rownames(mammTrait1) <- mammTrait$MammalName
colnames(mammTrait1)
# only species with trait data
mammTrait2 <- mammTrait1[complete.cases(mammTrait1),]
# number of species lost... 
dim(mammTrait1)-dim(mammTrait2)
# Make an RDA to scale the variance in diet
dietRDA <- vegan::rda(mammTrait2[,1:10])

dietMamm <- scores(dietRDA)$sites

mammTrait3 <- mammTrait1[complete.cases(mammTrait1),][c("Mass")]
mammTrait3 <- data.frame(mammTrait3,dietMamm)
mammTrait4 <- mammTrait3[!mammTrait3$Mass == -9999,]
dim(mammTrait3)-dim(mammTrait4) # species lost


mammTrait4$Mass <- log(mammTrait4$Mass )
colnames(mammTrait4) <- c( "Body_Mass", "Diet1", "Diet2")


# subsetting by interactions

N <- xtabs( ~ PALM + FRUGIVORE, 
            data = palmMetaNet[
              palmMetaNet$biogeographicRegion == "Neotropics" & 
                palmMetaNet$frugClass == "MAMMAL" ,])
dim(N) # network size (P x F)
table(N==1)
length(rownames(palmTraits1))
# subsetting for species found in the metanetwork
palmTraits2 <- na.omit(palmTraits1[match(rownames(N), rownames(palmTraits1)),])
names(palmTraits2)
# subset to the species present in the metanetwork 
mammTrait5 <- na.omit(mammTrait4[match(colnames(N),
                                       rownames(mammTrait4)),])

length(mammTrait4)

################################################################################
# setup data matrices for RLQ analysis based on the pooled binary metaweb
# (as in Albrecht et al. 2018, Nat. Comm. 9:3177)
################################################################################


pTRLQ <- palmTraits2
mTRLQ <- mammTrait5

nRLQ <- N[rownames(pTRLQ),rownames(mTRLQ)]
nRLQ <- decostand(nRLQ, "pa")

dim(nRLQ)
dim(pTRLQ)
dim(mTRLQ)



################################################################################
# perform correspondence and principal components analyses
################################################################################
coa <- dudi.coa(as.data.frame(nRLQ), scannf = FALSE, nf = 4)
duP <- dudi.pca(as.data.frame(pTRLQ), scannf = FALSE, 
                center = TRUE, 
                scale = TRUE, 
                nf = 4,
                row.w = coa$lw)
duA <- dudi.pca(as.data.frame(mTRLQ),
                scannf = FALSE, center = TRUE, 
                scale = TRUE, nf = 4, row.w = coa$cw)
################################################################################
# perform RLQ analysis
################################################################################
RLQ <- rlq(duP, coa, duA, scannf = FALSE, nf = 4)
plot(RLQ)
# Calculate contribution of axis to variance
RLQ$eig[1]/sum(RLQ$eig)

summary(RLQ)



################################################################################
# perform permutation tests
################################################################################
# Global test for association between trait spaces of plants and animals
################################################################################
trRLQ <- fourthcorner2(as.data.frame(pTRLQ),
                       as.data.frame(nRLQ),
                       as.data.frame(mTRLQ),
                       p.adjust.method.G = "fdr",
                       modeltype = 6, nrepet = 9999)$trRLQ

plot(trRLQ$plot$hist)
abline(v = trRLQ$obs,col =  "red")


trRLQ$pvalue
trRLQ$obs
trRLQ$expvar

################################################################################
# Test for axis-specific correlations between plant and animal trait spaces
# NOTE: only tests of 1st vs. 1st and 2nd vs. 2nd axis are meaningful
################################################################################
fcRLQ <- fourthcorner.rlq(RLQ, type = "axes",
                          p.adjust.method.G = "fdr",
                          p.adjust.method.D = "fdr",
                          p.adjust.D = "levels",
                          modeltype = 6, nrepet = 9999)
fcRLQ$tabD$obs[c(1,6)]
fcRLQ$tabD$expvar[c(1,6),]

plot(fcRLQ$tabG$plot[[6]]$hist)
abline(v=fcRLQ$tabG$obs[fcRLQ$tabG$pvalue<0.05], col = "red")



fcRLQ$tabG$names[which(fcRLQ$tabG$pvalue<0.06)]
################################################################################
# Test for association of animal traits with first and second axis of animal trait space
################################################################################
qAxes <-  fourthcorner.rlq(RLQ, type = "Q.axes",
                           p.adjust.method.G = "fdr",
                           p.adjust.method.D = "fdr",
                           p.adjust.D = "levels",
                           modeltype = 6, nrepet = 9999)

qAxes$tabD$obs[qAxes$tabD$pvalue < 0.1]

rAxes$tabD$obs[rAxes$tabD$pvalue < 0.06]


################################################################################
# Test for association of plant traits with first and second axis of plant trait space
################################################################################
rAxes <- fourthcorner.rlq(RLQ, type = "R.axes",
                          p.adjust.method.G = "fdr",
                          p.adjust.method.D = "fdr",
                          p.adjust.D = "levels",
                          modeltype = 6, nrepet = 9999)
rAxes
################################################################################


################################################################################
################################################################################
# the scores from the RLQ analysis can be used to calculate functional diversity
# indices on each elevation and conduct null model analyses
################################################################################
# extract normed species scores of plants and animals for further analysis
# these scores represent the position of plants and animals in the
# respective trait spaces
################################################################################
# Animals - extract normed species scores in trait space & create distance matrix
aDis <- dist(RLQ$mQ)
# Plants - extract normed species scores in trait space & create distance matrix
RLQ$mR
pDis <- dist(RLQ$mR)

