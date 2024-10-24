#setwd("/Users/jasminalbert/Documents/HALLETTLAB/mhav/code")
#datLoc <- "../data/"
seedsLoc <- paste0(datLoc, "competition_all_seeds_2021.csv")
seeds <- read.csv(seedsLoc, row.names=1, stringsAsFactors=T)
seeds$id <- 1:nrow(seeds)
sp <- c("AVBA", "TACA")
mas <- seeds[seeds$phyto%in%sp&seeds$background%in%c(sp,"Control"),]
#treat NAs in phyto data as zero
startpdat <- which(names(mas)=="tot_stems")
endpdat <- which(names(mas)=="output")
length(mas[,startpdat:endpdat][is.na(mas[,startpdat:endpdat])])==sum(is.na(mas[,startpdat:endpdat]))
mas[,startpdat:endpdat][is.na(mas[,startpdat:endpdat])] <- 0
sum(is.na(mas[,startpdat:endpdat]))==0

str(mas)
#saveRDS(mas,paste0(datLoc,"MHAVseeds.RDS"))
bg <- unique(mas$background)
acol <- "#b37746ff"
mcol <- "#80d850ff"
cols <- c(AVBA=acol,TACA=mcol)
comps <- expand.grid(fc=unique(mas$phyto),bg=unique(mas$background))
fmax <- max(mas$tot_seeds,na.rm=T)
bmax1 <- max(mas$mean_density_halfm2,na.rm=T)
bmax2 <- max(mas$max_density_halfm2,na.rm=T)
ylim=c(0,fmax)
x <- seq(0,bmax1,2)		
stats <- as.list(comps[,1]);names(stats)<-rownames(comps)
bgend <- 0
for (p in rownames(comps)){
	fc <- comps[p,"fc"]
	bg <- comps[p,"bg"]
	pair <- mas$phyto==fc& mas$background==bg
	fdat <- mas$tot_seeds[pair]
	bdat1 <- mas$mean_density_halfm2[pair]
	#bdat2 <- mas$max_density_halfm2[pair]
	main <- paste0("fc:",fc,"-bg:",bg)

	plot(bdat1,fdat, main=main, xlim=c(0,bmax1),ylim=ylim,col=cols[fc])
	#plot(bdat2,fdat, main=main, xlim=c(0,bmax2),ylim=ylim)
	if (bg!="Control"){
	model <- summary(lm(fdat~bdat1))
	m <- model$coefficients[2,1]
	b <- model$coefficients[1,1]
	pv <- model$coefficients[2,4]
	r <- model$r.squared
	stats[[p]] <- c(m,b,pv,r)
	y <- m*x+b
	lines(x,y)
	text(0,fmax*.83,paste0("y=",signif(m,3),"x+", signif(b,3),"\np=", signif(pv,3),", R2=",signif(r,3)), adj=c(0,0))
	if (bg!=bgend){
		plot(fdat,bdat1, main=main, xlim=c(0,bmax1),ylim=ylim,col=cols[fc])
	} else if(bg==bgend){points(fdat,bdat1,col=cols[fc])}
	bgend <- bg
	}
}



