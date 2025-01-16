getdata <- function(data, sp, nut="all", ppt="all"){
	sps <- as.character(unique(data$phyto))
	dat <- data[data$phyto==sp,]
	if(nut!="all"){
		dat <- dat[dat$nut_trt==nut,]
	}
	if(ppt!="all"){
		dat <- dat[dat$ppt_trt==ppt,]
	}
	na.rows <- rowSums(apply(dat,2,is.na))>0
	dat <- dat[!na.rows,]
	
	#create model variables
	## set population context for each species as seeds in
	av <- as.integer(dat$AVBAIn)
	mh <- as.integer(dat$TACAIn)
	
	## set fecundity as seeds out of focal sp
	spSeeds <- which(substring(names(dat),1,4)==sp)
	names(spSeeds) <- substring(names(dat)[spSeeds],5)
	Fecundity <- as.integer(dat[,spSeeds["Out"]]) 
	
	## set intraspecific species
	intra <- as.integer(dat[,spSeeds["In"]])
	Plot <- dat$block
	
	## number of observations
	N <- length(Fecundity)
	P <- length(unique(Plot))
	
	## germination rates
	ag <- 0.92
	mg <- 0.87
	
	res <- list("N"=N,"Fecundity"=Fecundity, "intra"=intra,"av"=av,"mh"=mh, "P"=P,"Plot"=Plot,"mg"=mg,"ag"=ag)
	list2env(res,.GlobalEnv)
	return(res)
}
#rm("N","Fecundity","intra","av","mh","P","Plot","mg","ag")
#getdata(data,"TACA",nut="XC",ppt="XC")
#ls()








