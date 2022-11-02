# Modelling FTA in relation to climate 
# Gabriel Muñoz
# Nov. 22


# recreate variables 

FRA_palms <- ZscoresINB
FRA_mammals <-ZscoresINBm[match(names(ZscoresINB), names(ZscoresINBm))]

FEA_palms <- ZscoresIND
FEA_mammals <- ZscoresINDm[match(names(ZscoresIND), names(ZscoresINDm))]


# visualize asymmetry

par(mfrow = c(1,2))
plot(FRA_palms~FRA_mammals, 
     xlim = c(0,8),
     ylim = c(0,8), 
     frame = F)
abline(a=0,b=1, lwd = 2)
abline(h=0,v=0, lty= 2)

plot(FEA_palms~FEA_mammals,
     xlim = c(-2,5),
     ylim = c(-2,5), 
     frame = F)
abline(a=0,b=1, lwd =2)
abline(h=0,v=0, lty= 2)

# download climatic data
WCLim <- raster::getData("worldclim", var="bio",res=10)
# function to cut to the neotropics
cropMask <- function(raster,prov){
  ## crop and mask
  r2 <- crop(raster, extent(prov))
  r3 <- mask(r2, prov)
  return(r3)
}


WCLim <- cropMask(WCLim, neotropic)

# Separate variables of interest
Temp <- WCLim[[1]]
Prec <- WCLim[[12]]
PrecSe <- WCLim[[15]]
IsoTer<- WCLim[[3]]
TempSeaso<- WCLim[[14]]

######
# Model climatic data 
####
Temp <- aggregate(Temp, 1/0.17)
Prec <- aggregate(Prec, 1/0.17)
PrecSe <- aggregate(PrecSe, 1/0.17)
TempSeaso <- aggregate(TempSeaso, 1/0.17)


Assemblages <- SpatialPoints(na.omit(data.frame(apply(stringr::str_split(names(FRA_mammals),
                                                                         "_", simplify = T), 
                                                      2, as.numeric))), 
                             proj4string = ecoRegi@proj4string)

TempAss <- extract(Temp,Assemblages)
PrecAss <- extract(Prec,Assemblages)
TempSeAss <- extract(TempSeaso,Assemblages)
PrecSeAss <- extract(PrecSe,Assemblages)

#####
# Make it all into a single model
####

# spatial gradients  
climatt <- data.frame("Te" = TempAss,
                      "P" = PrecAss, 
                      "TS" = TempSeAss,
                      "PS" = PrecSeAss)
climatt <- apply(climatt, 2, decostand, "standardize")
climatt <- data.frame(climatt)
climatt$Dom <- extract(neotropic,Assemblages)$Dominions

## Fit linear models
toRem <- which(is.na(names(FRA_mammals)))


# Fit linear mixed models 
# Mammal functional richness
M_MamSize_clim <- lm(scale(FRA_mammals[-toRem])~ Te+P+TS+PS , data = climatt)
M_MamSize_dom <- lm(scale(FRA_mammals[-toRem])~ Dom , data = climatt)
M_MamSize_add <- lm(scale(FRA_mammals[-toRem])~ Te+P+TS+PS+Dom , data = climatt)
M_MamSize_mix <- lm(scale(FRA_mammals[-toRem])~ Te*Dom+P*Dom , data = climatt)

AIC <- AIC(M_MamSize_clim,M_MamSize_dom,M_MamSize_add,M_MamSize_mix)
AIC[order(AIC$AIC),]
AICcmodavg::aictab(list(M_MamSize_clim,
                        M_MamSize_dom,M_MamSize_add,M_MamSize_mix))

RsquareAdj(M_MamSize_clim)
RsquareAdj(M_MamSize_dom)
RsquareAdj(M_MamSize_add)
RsquareAdj(M_MamSize_mix)
anova(M_MamSize_mix)

# Palms functional richness

P_PalmSize_clim <- lm(scale(FRA_palms[-toRem])~ Te+P+TS+PS , data = climatt)
P_PalmSize_dom <- lm(scale(FRA_palms[-toRem])~ Dom , data = climatt)
P_PalmSize_add <- lm(scale(FRA_palms[-toRem])~ Te+P+TS+PS+Dom , data = climatt)
P_PalmSize_mix <- lm(scale(FRA_palms[-toRem])~ Te*Dom+P*Dom , data = climatt)

AICcmodavg::aictab(list(P_PalmSize_clim,P_PalmSize_dom,P_PalmSize_add,P_PalmSize_mix))

RsquareAdj(P_PalmSize_clim)
RsquareAdj(P_PalmSize_dom)
RsquareAdj(P_PalmSize_add)
RsquareAdj(P_PalmSize_mix)
summary(P_PalmSize_mix)


## Mammal functional evenness


M_MamEven_clim <- lm(scale(FEA_mammals[-toRem])~ Te+P+TS+PS , data = climatt)
M_MamEven_dom <- lm(scale(FEA_mammals[-toRem])~ Dom , data = climatt)
M_MamEven_add <- lm(scale(FEA_mammals[-toRem])~ Te+P+TS+PS+Dom , data = climatt)
M_MamEven_mix <- lm(scale(FEA_mammals[-toRem])~ Te*Dom+P*Dom , data = climatt)

AICcmodavg::aictab(list(M_MamEven_clim,M_MamEven_dom,M_MamEven_add,M_MamEven_mix))

vegan::RsquareAdj(M_MamEven_clim)
RsquareAdj(M_MamEven_dom)
RsquareAdj(M_MamEven_add)
RsquareAdj(M_MamEven_mix)

# Palm functional evenness

P_PalmEven_clim <- lm(scale(FEA_palms[-toRem])~ Te+P+TS+PS , data = climatt)
P_PalmEven_dom <- lm(scale(FEA_palms[-toRem])~ Dom , data = climatt)
P_PalmEven_add <- lm(scale(FEA_palms[-toRem])~ Te+P+TS+PS+Dom , data = climatt)
P_PalmEven_mix <- lm(scale(FEA_palms[-toRem])~ Te*Dom+P*Dom , data = climatt)

AICcmodavg::aictab(list(P_PalmEven_clim,P_PalmEven_dom,P_PalmEven_add,P_PalmEven_mix))
easystats::model_dashboard(P_PalmEven_mix)

RsquareAdj(P_PalmEven_clim)
RsquareAdj(P_PalmEven_dom)
RsquareAdj(P_PalmEven_add)
RsquareAdj(P_PalmEven_mix)

# Visualize results 
sjPlot::plot_models(M_MamSize_mix,P_PalmSize_mix,
                    M_MamEven_mix, P_PalmEven_mix,
                    grid = T, p.shape = T)

## Calculate FTA and model its variation with climate


AsyClimate <- data.frame("SizeAsy"= c(FRA_mammals[-toRem]-FRA_palms[-toRem]),
                         "EveAsy" = c(FEA_mammals[-toRem]-FEA_palms[-toRem]),
                         "H2" = H2_mean$value[match(names(ZscoresIND),
                                                    H2_mean$L1)][-toRem],
                         climatt)


summary(lm(AsyClimate$SizeAsy~ 
          AsyClimate$EveAsy))
   
plot(AsyClimate$EveAsy~ 
     AsyClimate$Te, 
     pch = 16, 
     col = f(AsyClimate$Te,9, "Reds", T))

abline(h=0,v=0, lty = 2)
dev.off()

SizeAsy_clim <- lm(SizeAsy~ Te+P+TS+PS , data = AsyClimate)
SizeAsy_clim_dom <- lm(SizeAsy~ Dom , data = AsyClimate)
SizeAsy_clim_add <- lm(SizeAsy~ Te+P+TS+PS+Dom , data = AsyClimate)
SizeAsy_clim_mix <- lm(SizeAsy~ Te*Dom+P*Dom , data = AsyClimate)


AICcmodavg::aictab(list(SizeAsy_clim,SizeAsy_clim_dom,SizeAsy_clim_add,SizeAsy_clim_mix))

RsquareAdj(SizeAsy_clim)
RsquareAdj(SizeAsy_clim_dom)
RsquareAdj(SizeAsy_clim_add)
RsquareAdj(SizeAsy_clim_mix)

PEvenAsy_clim <- lm(EveAsy~ P*Dom , data = AsyClimate)
EvenAsy_clim <- lm(EveAsy~ Te+P+TS+PS , data = AsyClimate)
EvenAsy_clim_dom <- lm(EveAsy~ Dom , data = AsyClimate)
EvenAsy_clim_add <- lm(EveAsy~ Te+P+TS+PS+Dom , data = AsyClimate)
EvenAsy_clim_mix <- lm(EveAsy~ Te*Dom+P*Dom , data = AsyClimate)

AICcmodavg::aictab(list(EvenAsy_clim,EvenAsy_clim_dom,EvenAsy_clim_add,EvenAsy_clim_mix))

RsquareAdj(EvenAsy_clim)
RsquareAdj(EvenAsy_clim_dom)
RsquareAdj(EvenAsy_clim_add)
RsquareAdj(EvenAsy_clim_mix)


sjPlot::plot_models(SizeAsy_clim_mix, 
                    EvenAsy_clim_mix, 
                    grid = T, p.shape = T)
