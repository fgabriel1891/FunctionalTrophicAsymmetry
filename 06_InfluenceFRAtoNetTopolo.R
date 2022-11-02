### Relating functional asymmetry to Network Topology


toRDA <- AsyClimate[complete.cases(AsyClimate),]
head(toRDA)


leg <- unique(data.frame("Dom" = as.factor(neotropic$Dominions),
       "col" = f(as.numeric(as.factor(neotropic$Dominions)), 7, "Set2")))

leg <- leg[!is.na(leg$col),]
par(mar = c(0,0,0,0))
plot(neotropic, 
     lwd = 0.3, col = f(as.numeric(as.factor(neotropic$Dominions)), 7, "Set2"))

legend("bottomleft", 
       pch = 22,
       cex = 1.5,
       title = "Biogeographic region",
       bty = "n",
       legend = leg$Dom, 
       pt.bg = leg$col)

# fit a multiple regression 
MR <- lapply(1:7, function(i) 
  lm(H2~abs(SizeAsy) + abs(EveAsy) ,toRDA[toRDA$Dom == leg$Dom[i],]))

for(i in 1:7){
  car::avPlots(MR[[i]], 
               id = F,
               col = leg$col[i],
               pch = 16,
               col.lines = "black",
               lwd = 3,
               xlim = c(-2,2),
               ylim = c(-0.2,0.2),
               marginal.scale=F,
               main = leg$Dom[i])
  }


