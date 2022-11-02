### Calculating palm animal cooccurrence at 1degree resolution 
## Gabriel Muñoz 
## Nov 2022


## Data on species traits

## Palm-frugivore data
intDataRaw <- read.csv("DATA/PalmDryadRepo/PalmFrugDatasetOCT2018.csv")
## Palm trait data 
palmTraits <- read.csv("DATA/PalmTraits_1.0.txt", sep = "\t", stringsAsFactors = F)
# conserve relevant columns
palmTraitsRel <- palmTraits[names(palmTraits)[c(1,7,8,12,13,14,19,20,21,22,23,26)]]






## Animal trait data 
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


birdsTr <- read.table("DATA/VAR/BIRDS/BirdFuncDat.txt", sep = )
# coercing bird and mammal datasets 
animalTrait <-  rbind(
  mammElt[names(birdsTr)[names(birdsTr)%in% names(mammElt)]],
  birdsTr[names(mammElt)[names(mammElt)%in% names(birdsTr)]]
)

names(animalTrait)
####
names(palmTraits)[1] <- "PALM"
names(animalTrait)[1] <- "FRUGIVORE"




library(raster)
library(rgdal)
## load palm distributional data 
# get shapefilenames 
filenames <- list.files("DATA/VAR/Palms/Shapefiles/")[grep(".shp", list.files("DATA/VAR/Palms/Shapefiles/"))]
# prune .shp.xml filetypes
filenames <- filenames[-grep(".xml", filenames)]
# recreate path
filepath <- paste0("DATA/VAR/Palms/Shapefiles/", filenames)

# remove shapes with errors
filenames <-filenames[-grep("Geonoma_arundinacea.shp", filenames)]
filenames <-filenames[-grep("TEMPLATE_WGS84.shp", filenames)]
filepath <- paste0("DATA/VAR/Palms/Shapefiles/", filenames)

## batch load the shapes
PalmPolygon <- list()

for(i in 1:length(filepath)){
  print(filepath[i])
  print(i)
  PalmPolygon[[i]] <- raster::shapefile(filepath[i])
  cat("\014")
  
}


library(raster)
# batch rasterize the polygons

palm.r <- list()
r <- list()
for(i in 1:length(filepath)){
  print(filepath[i])
  print(i)
  r[[i]] <- raster()
  extent(r[[i]]) <- extent(PalmPolygon[[i]])
  palm.r[[i]] <- rasterize(PalmPolygon[[i]], r[[i]])
  cat("\014")
}

## set a similar resolution for all layers
# make the reference layer
f <- raster(ncol = 140, nrow = 240)
extent(f) <- c(-120,-30,-60,60)



# for(i in 1:length(filepath)){
#   print(filepath[i])
#   print(i)
#   palm.r[[i]] <- resample(palm.r[[i]], f, "bilinear")
#   cat("\014")
# }
# 
# palm.stack <- stack(palm.r)
# 
# 
library(stringr)
spNames <- stringr::str_split(filenames, pattern = ".shp", simplify = T)[,1]
spNames <- paste(str_split(spNames, "_",  simplify = T)[,1],
                 str_split(spNames, "_",  simplify = T)[,2])
comFile <- c()






## iterate over all the shapefiles
PalmPresenceDat <- lapply(1:length(palm.r),
                          function(x) getPalmValues(palm.r[[x]], 
                                                    spNames[x]))


## bind the nested data.frames into a single one. 
PalmPresenceDat <- plyr::rbind.fill(PalmPresenceDat)
## Round to 1x1 degree resolution 
PalmPresenceDat$x.R <- round(PalmPresenceDat$x)
PalmPresenceDat$y.R <- round(PalmPresenceDat$y)
## Make site unique id
PalmPresenceDat$unID <- paste0(PalmPresenceDat$x.R, "_",PalmPresenceDat$y.R)



# Save image! 

# saveRDS(PalmPresenceDat, "Scripts/FunctEcology/PalmPointDat.RDS")

#######
## Important! check names of palms--- resolve later
#######

palmNamesNonMatch <- unique(PalmPresenceDat$spName[!PalmPresenceDat$spName
                                                   %in% palmTraitsRel$SpecName])
length(palmNamesNonMatch) # 52 species with range info are not present in trait dataset

######
##  Animal interaction niche.  
#######

# build a new dataset with the 1x1 distribution and functional diversity of traits 

TraitDistPalm <- cbind(PalmPresenceDat,palmTraitsRel[match(PalmPresenceDat$spName,
                                                           palmTraitsRel$SpecName),])

## prune the dataset to include only traits for species with matching distribution data 
palmTraitPrun <- palmTraitsRel[palmTraitsRel$SpecName %in% PalmPresenceDat$spName,]
PalmPresenceDatPrun <- PalmPresenceDat[PalmPresenceDat$spName %in% palmTraitPrun$SpecName,]
PalmPresenceDatPrun <- droplevels(PalmPresenceDatPrun)

#########
# Plant interaction niche
########

## select species with some degree of frugivory
frugTraits <- animalTrait[animalTrait$Diet.Fruit > 0,]
frugTraits <- frugTraits[!c(names(frugTraits) %in% c("Diet.Source","Diet.Certainty",
                                                     "BodyMass.Source", "BodyMass.SpecLevel"))]


## load species distribution data

mammShape <- shapefile("AnimalDistribution/Mammal/TERRESTRIAL_MAMMALS.shp")
# get the shapes corresponding with animal frugivores with we have trait information data 
mammShape <- mammShape[mammShape@data$binomial %in% frugTraits$FRUGIVORE,]
# subset those shapes so it only contains america 
# get the lower limit (lat)
Corxmax <- sapply(1:length(mammShape$binomial), 
                  function(x) extent(bbox(subset(mammShape, mammShape$binomial == mammShape$binomial[[x]])))@xmax)
# subset only for american species
AmMammal <- lapply(which(Corxmax <0), function(x) subset(mammShape, mammShape$binomial == mammShape$binomial[[x]]))


# batch rasterize the polygons

amMammal.r <- list()
r <- list()
for(i in 1:length(AmMammal)){
  print(i)
  r[[i]] <- raster()
  extent(r[[i]]) <- extent(AmMammal[[i]])
  amMammal.r[[i]] <- rasterize(AmMammal[[i]], r[[i]])
  cat("\014")
}

## set a similar resolution for all layers
# make the reference layer
f2 <- raster(ncol = 140, nrow = 240)
extent(f2) <- c(-120,-30,-60,60)

# get the species names from the shapefiles
MammSp <- sapply(1:length(AmMammal), function(x) AmMammal[[x]]@data$binomial[1])

comFile <- c()


## iterate over all the shapefiles
MammPresenceDat <- lapply(1:length(amMammal.r), 
                          function(x) getPalmValues(amMammal.r[[x]], MammSp[x]))

## bind the nested data.frames into a single one. 
MammPresenceDat <- plyr::rbind.fill(MammPresenceDat)
## Round to 1x1 degree resolution 
MammPresenceDat$x.R <- round(MammPresenceDat$x)
MammPresenceDat$y.R <- round(MammPresenceDat$y)
## Make site unique id
MammPresenceDat$unID <- paste0(MammPresenceDat$x.R, "_",MammPresenceDat$y.R)

length(unique(MammPresenceDat$spName))
# save image

# saveRDS(MammPresenceDat, "Scripts/FunctEcology/MammalPointDat.RDS")

#######
## Important! check names of mammal species--- resolve later all included!
#######

MammalNamesNonMatch <- unique(MammPresenceDat$spName[!MammPresenceDat$spName %in% frugTraits$FRUGIVORE])

######
##  Palm interaction niche.  
#######


# build a new dataset with the 1x1 distribution and functional diversity of traits 

TraitDistMammal <- cbind(MammPresenceDat,frugTraits[match(MammPresenceDat$spName, frugTraits$FRUGIVORE),])

## prune the dataset to include only traits for species with matching distribution data (ideally this step should be removed, by matching all names properly... later on...)
MammalTraitPrun <- frugTraits[frugTraits$FRUGIVORE %in% MammPresenceDat$spName,]
## remove those species for we don't have trait data!
# (ideally this step should be removed, 
# by matching all names properly and gathering more trait info. later on...)

MammPresenceDatPrun <- MammPresenceDat[MammPresenceDat$spName %in% frugTraits$FRUGIVORE,]
MammPresenceDatPrun <- droplevels(MammPresenceDatPrun)







