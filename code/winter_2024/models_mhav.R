datLoc <- "../data/"
ATdat <- readRDS(paste0(datLoc, "dataAvbaTaca.RDS"))

names(ATdat)[which(names(ATdat)=="max_density_halfm2")] <- "bgSeedIn"
names(ATdat)[which(names(ATdat)=="mean_density_halfm2")] <- "bgSeedOut"
names(ATdat)[which(names(ATdat)=="tot_seeds")] <- "fcSeedOut"
ATdat$fcSeedIn <- 3

ID <- paste(ATdat$block, ATdat$plot, ATdat$subplot, ATdat$phyto, sep="_")
ATdat$ID <- ID; rownames(ATdat) <- ID
subID <- paste(ATdat$block, ATdat$plot, ATdat$subplot, sep="_"); subID <- unique(subID)

#need to reshape data with TACA seedin/out and AVBA seedin/out (regardless of role)

Tphyto <- ATdat[ATdat$phyto=="TACA",]
Aphyto <- ATdat[ATdat$phyto=="AVBA",]
dim(ATdat);dim(rbind(Tphyto, Aphyto))
rownames(Aphyto) <- rownames(Tphyto) <- subID
ATlist <- list(Aphyto, Tphyto)
phytos <- c("av","tc")
names(ATlist) <- phytos

datnames <- c( paste0(phytos[1],c("SeedIn", "SeedOut")), paste0(phytos[2],c("SeedIn", "SeedOut")))
#datnames <- c( "SeedIn", "SeedOut")
tempmat <- matrix(0, nrow=length(subID), ncol=length(datnames), dimnames=list(subID, datnames))
outdf <- data.frame(tempmat)
outlist <- list(outdf, outdf); names(outlist) <- phytos 

for (s in phytos){
	seedMets <- paste0(s,c("SeedIn", "SeedOut"))
	seedMetsBg <- paste0(phytos[phytos!=s],c("SeedIn", "SeedOut"))
	for (i in subID){
		sample <- ATlist[[s]][i,]
		outlist[[s]][i,seedMets] <- sample[,c("fcSeedIn", "fcSeedOut")] 
		
		if (sample$background == sample$phyto){
			outlist[[s]][i,seedMets] <- outlist[[s]][i,seedMets] + sample[,c('bgSeedIn', 'bgSeedOut')]
		} else if (sample$background == "Control"){
			outlist[[s]][i,seedMetsBg] <- c(0,0) 
		} else {
			outlist[[s]][i,seedMetsBg] <- sample[,c("bgSeedIn", "bgSeedOut")]}
	}	
	rownames(outlist[[s]]) <- paste(rownames(outlist[[s]]), s, sep="_")
}
ATseeds <- rbind(outlist[[1]],outlist[[2]])




ag <- .9
tg <- .9

mT <- as.formula(log(tcSeedOut +1) ~  log(tg*(tcSeedIn+1)*exp(log(lambda)-log((1+aiA*(avSeedIn+1)*ag+aiT*(tcSeedIn+1)*tg)))))

mA <- as.formula(log(avSeedOut +1) ~  log(ag*(avSeedIn+1)*exp(log(lambda)-log((1+aiT*(tcSeedIn+1)*tg+aiA*(avSeedIn+1)*ag)))))

#fit model
AVoutput <- as.data.frame(matrix(nrow = 0, ncol = 7))
names(AVoutput) = c("estimate", "se", "t", "p", "params", "treatment", "species")
m1out <- nlsLM(mA, start=list(lambda=1, aiT = .01, aiA=.01),control=nls.lm.control(maxiter=500),lower = c(0, 0, 0), upper = c(200, 1, 1), trace=T, data = ATseeds)

outreport <- as.data.frame(summary(m1out)$coef[1:3, 1:4])
names(outreport) = c("estimate", "se", "t", "p")
outreport$params <- row.names(outreport)
outreport$treatment <- treatments[i]
outreport$species <- "Avena"
AVoutput <- rbind(AVoutput, outreport)

TCoutput <- as.data.frame(matrix(nrow = 0, ncol = 7))
names(TCoutput) = c("estimate", "se", "t", "p", "params", "treatment", "species")
mTout <- nlsLM(mT, start=list(lambda=1, aiT = .01, aiA=.01),control=nls.lm.control(maxiter=500),lower = c(0, 0, 0), upper = c(200, 1, 1), trace=T, data = ATseeds)

outreport <- as.data.frame(summary(m1out)$coef[1:3, 1:4])
names(outreport) = c("estimate", "se", "t", "p")
outreport$params <- row.names(outreport)
outreport$treatment <- treatments[i]
outreport$species <- "Taca"
TCoutput <- rbind(TCoutput, outreport)






