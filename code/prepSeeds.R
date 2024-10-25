#function to seeds in seeds out data
stem2seedBg <- function(mas){
	mhOutputIntra <- mas[mas$phyto=="TACA" & 	mas$role=="competitor"&mas$output>0,]$output
	avOutputIntra <- mas$output[mas$phyto=="AVBA" & mas$role=="competitor"&mas$output>0]
	mhmn<-mean(mhOutputIntra);mhsd<-sd(mhOutputIntra)
	avmn<-mean(avOutputIntra);avsd<-sd(avOutputIntra)
	avmn <- 25;sd <- 25
	avran <- exp(rnorm(500,log(avmn),log(avsd)*.3))
	avran <- avran[avran<65&avran>15]
	mhran <- rnorm(500,mhmn,mhsd*.5)
	rvarstem2seeds <- cbind(TACA=mhran[1:length(avran)],AVBA=avran)
	return(rvarstem2seeds)
}

prepSeeds <- function(seeddat){
	#stem to seeds distribution for estimating bg seeds out
	rvarstem2seeds <- stem2seedBg(seeddat)
	
	sp <- as.character(unique(seeddat$phyto))
	seedstore <- data.frame(matrix(ncol=2,dimnames= list(NULL,sp)))
	seedres <- list(In=seedstore,Out=seedstore)
	
	for (i in 1:nrow(seeddat)){
		sdat <- seeddat[i,]
		isp <- as.character(sdat$phyto)
		seedsin <- c(phy=3,bg=sdat$max_density_halfm2/6)
		sam <- sample(1:nrow(rvarstem2seeds),1)
		seedsout <- c(phy=sdat$tot_seeds,bg=(sdat$mean_density_halfm2/6)* rvarstem2seeds[[sam,isp]])
		if(sdat$background==isp){
			phyto <- sum(seedsin)
			background <- 0
			phyto_o <-  sum(seedsout)
			background_o <- background
		} else {
			phyto <- seedsin["phy"]; 
			background <- seedsin["bg"]
			phyto_o <- seedsout['phy']
			background_o <- seedsout['bg']
		}
		seedres$In[i,isp] <- phyto
		seedres$In[i,sp!=isp] <- background 	
		seedres$Out[i,isp] <- phyto_o
		seedres$Out[i,sp!=isp] <- background_o	
	}
	seed <- cbind(seedres$In,seedres$Out)
	names(seed) <- paste0(names(seed),rep(names(seedres),each=2))
	return(seed)
}
#seed2 <- prepSeeds(seeddat)
#seed2==seed





