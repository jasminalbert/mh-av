#working directing is /mhav/code
datLoc <- "../data/"
seedsLoc <- paste0(datLoc, "competition_all_seeds_2021.csv")
seeds <- read.csv(seedsLoc, row.names=1, stringsAsFactors=T)
head(seeds)

## notes from ashely ##
#output = seeds/stems
#mean_density_halfm2 = observed density of background competitor in the field (stem counts). some missing
#max_density_halfm2 = maximum possible density calculated based on seed mass (seeding at 8g/m2)

unique(seeds$nut_trt)
unique(seeds$ppt_trt)
unique(seeds$phyto)
unique(seeds$background)
unique(seeds$role)
tacabg <- seeds[seeds$background=="TACA",]
nrow(tacabg);head(tacabg)
par(mar=c(4,4,1,3))
plot(seeds$phyto, seeds$role)
plot(tacabg$phyto, tacabg$role)

plot(tacabg$output, tacabg$mean_density_halfm2, col=tacabg$nut_trt);legend("topright", legend=levels(tacabg$nut_trt), fill=factor(levels(tacabg$nut_trt)))
plot(tacabg$max_density_halfm2, tacabg$mean_density_halfm2)
plot(tacabg$output, tacabg$max_density_halfm2, col=tacabg$phyto); legend("topright", legend=levels(tacabg$phyto), fill=tacabg$phyto)
plot(tacabg$output, tacabg$mean_density_halfm2, col=tacabg$phyto);legend("topright", legend=levels(tacabg$phyto), fill=factor(levels(tacabg$phyto)))
plot(tacabg$output, tacabg$mean_density_halfm2, col=tacabg$ppt_trt);legend("topright", legend=levels(tacabg$ppt_trt), fill=factor(levels(tacabg$ppt_trt)))



#only TACA and AVBA
target <- c("AVBA", "TACA", "Control")
atseeds <- seeds[seeds$phyto%in%target & seeds$background%in%target,]
nrow(atseeds);head(atseeds)
atcntrl <- atseeds[atseeds$nut_trt=="N" & atseeds$ppt_trt=="XC",]

#plots
pdf("../figures/dataPointsN.pdf")
par(mfrow=c(3,2), mgp=c(1.3,0.4,0), tcl=-.3, mar=c(2,3,0.5,1), oma=c(1,0,0,0), xpd=NA)
for (p in unique(atcntrl$phyto)){
	for (b in unique(atcntrl$background)){
		out <- atcntrl[atcntrl$phyto==p & atcntrl$background==b,]
		plot(out$output, out$mean_density_halfm2, xlab=paste(p,"output"), ylab=paste(b,"density"), bty='l', xlim=c(0,60), ylim=c(0,250))
		title(main=paste(p,"in",b), font.main=1,line=-1)
	}
}
dev.off()

saveRDS(atcntrl, paste0(datLoc,"dataAvbaTaca_N.RDS"))














