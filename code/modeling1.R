#setwd("/Users/jasminalbert/Documents/HALLETTLAB/mh-av/code")
library(minpack.lm)
datLoc <- "../data/"
filename <- "MHAVseeds.RDS"
mas <- readRDS(paste0(datLoc,filename))
#mas <- mas[mas$nut_trt!="N",] 	
	#if using max stems out, getting Ntrt is fine
str(mas)
max(mas$tot_stems)
#seedsIn =3, but max stems out is 4
#stems in stems out model or seedsind seedsout
#seeds in: have for phytos but background is in grams seeds
#seeds out:have for phytos but background is stems out
#stems in: no
#stems out: yes for both
#can i do seeds in seeds out but alpha are not symmetric units?
#so alpha_i,j is effects of stems of j - no seed in timet
#seeds in for phytos and bg to seeds out for phyto 
#then i wont even use the background data?

#seeds in for bg = 8g/m2
#subplots are 0.5x0.5 m2 = 0.25m2
#=8*.25 = 2g/m2
#MH one seed = 0.0058g (from JS)
2/0.0058 #=344.8276 seeds in
#AV: .015g
2/0.015 #=133.3333 seedsin
#AS expected 530.504 and 116.959 max stems out from seeding rate
#so that means she has one seed = 
2/530.504 #=0.00377 MH g/seed
2/116.959 #=0.01710001 AV g/seed  
#this is for the whole subplot but the plot was divided into six subsubplots?
#divide maxstemsout by 6 to get bg seeds in?
head(mas)
reshape(mas,direction="long")
mhOutputIntra <- mas[mas$phyto=="TACA" & mas$role=="competitor"&mas$output>0,]$output
avOutputIntra <- mas$output[mas$phyto=="AVBA" & mas$role=="competitor"&mas$output>0]
mhmn<-mean(mhOutputIntra);mhsd<-sd(mhOutputIntra)
avmn<-mean(avOutputIntra);avsd<-sd(avOutputIntra)
avmn <- 25
sd <- 25
hist(mhOutputIntra)
hist(avOutputIntra,10)
range(avOutputIntra)
hist(rnorm(500,mhmn,mhsd*.5))
range(rnorm(500,mhmn,mhsd*.5))
avran <- rnorm(500,avmn,avsd*.3)
avran <- exp(rnorm(500,log(avmn),log(avsd)*.3))
avran <- avran[avran<65&avran>15]
hist(avran,breaks=10)
rvarstem2seeds <- cbind(TACA=mhran[1:length(avran)],AVBA=avran)
head(rvarstem2seeds)
idend <- which(names(mas)=="subplot")
datid <- mas[,1:idend]
seeddat <- mas[,-(1:idend)];head(seeddat)
seeddat <- seeddat[,!names(seeddat)%in%c("phytonum", "role")]
dat <- data.frame(matrix(ncol=4,dimnames= list(NULL,c("mhIn","mhOut","avIn","avOut"))))
sp <- c("TACA","AVBA")
seedstore <- data.frame(matrix(ncol=2,dimnames= list(NULL,sp)))
seedres <- list(In=seedstore,Out=seedstore)
for (i in 1:nrow(seeddat)){
	sdat <- seeddat[i,]
	isp <- as.character(sdat$phyto)
	seedsin <- c(phy=3,bg=sdat$max_density_halfm2)
	sam <- sample(1:nrow(rvarstem2seeds),1)
	seedsout <- c(phy=sdat$tot_seeds,bg=sdat$mean_density_halfm2* rvarstem2seeds[[sam,isp]])
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
n<-sample(1:216,1);cbind(seedres$In,seedres$Out,seeddat)[n,]


seed <- cbind(seedres$In,seedres$Out)
names(seed) <- paste0(names(seed),rep(names(seedres),each=2))

ag <- .9
#m1E <- as.formula(log(ERseedout +1) ~  log(eg*(ERseedin+1)*exp(log(lambda)-log((1+aiE*(ERseedin+1)*eg+aiA*(AVseedin+1)*ag)))))

m1m <- as.formula(log(TACAOut +1) ~  log((TACAIn+1)*exp(log(lambda)-log(1+aiM*(TACAIn+1)+aiA*(AVBAIn+1)))))
MHoutput <- as.data.frame(matrix(nrow = 0, ncol = 6))
names(MHoutput) = c("estimate", "se", "t", "p", "params", "species")

m1out <- nlsLM(m1m, start=list(lambda=1, aiM = .01, aiA=.01),
                 lower = c(0, 0, 0), upper = c(200, 1, 1),
                 control=nls.lm.control(maxiter=500), trace=T,
                 data = seed[!is.na(seed$TACAOut),])
  
 outreport <- as.data.frame(summary(m1out)$coef[1:3, 1:4])
  names(outreport) = c("estimate", "se", "t", "p")
  outreport$params <- row.names(outreport)
  outreport$species <- "MH"
  MHoutput <- rbind(MHoutput, outreport)










