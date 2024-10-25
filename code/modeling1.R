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
#head(mas)
#reshape(mas,direction="long")

idend <- which(names(mas)=="subplot")
datid <- mas[,1:idend]
seeddat <- mas[,-(1:idend)];head(seeddat)
seeddat <- seeddat[,!names(seeddat)%in%c("phytonum", "role")]
dat <- data.frame(matrix(ncol=4,dimnames= list(NULL,c("mhIn","mhOut","avIn","avOut"))))

seeddat <- mas
#n<-sample(1:216,1);cbind(seedres$In,seedres$Out,seeddat)[n,]
seed <- prepSeeds(seeddat)
saveRDS(seed,"si-so.RDS")

ag <- .9
#m1E <- as.formula(log(ERseedout +1) ~  log(eg*(ERseedin+1)*exp(log(lambda)-log((1+aiE*(ERseedin+1)*eg+aiA*(AVseedin+1)*ag)))))

m1m <- as.formula(log(TACAOut +1) ~  log((TACAIn+1)*exp(log(lambda)-log(1+aiM*(TACAIn+1)+aiA*(AVBAIn+1)))))

#m1A <- as.formula(log(AVseedot +1) ~  log(ag*(AVseedin+1)*exp(log(lambda)-log((1+aiE*(ERseedin+1)*eg+aiA*(AVseedin+1)*ag)))))
m1A <- as.formula(log(AVBAOut +1) ~  log((AVBAIn+1)*exp(log(lambda)-log(1+aiM*(TACAIn+1)+aiA*(AVBAIn+1)))))
models <- list(mh=m1m,av=m1A)
modelOut <- fitModels(models,seed1)
modelOut==mo




plot(seed$TACAIn[seed$AVBAIn==0],seed$TACAOut[seed$AVBAIn==0])
plot(seed$TACAIn,seed$TACAOut)
plot((log(seed$TACAIn+1)+1)/(log(seed$AVBAIn+1)+1),log(seed$TACAOut+1))


sz <- round(log(seed$TACAOut+1))*.5+0.5
cl <- round(log(seed$TACAOut+1))+1
cols <- hcl.colors(20,"Viridis", alpha=0.1,T)[1:9]
plot(log(seed$TACAIn+1),log(seed$AVBAIn+1), cex=sz, col=cols[cl],pch=19, xlim=c(0,5), ylim=c(-0.5,4))
sz <- round(log(seed$AVBAOut+1))*.5+0.5
cl <- round(log(seed$AVBAOut+1))+1
cols <- hcl.colors(20,"Viridis", alpha=0.1,T)[1:7]
plot(log(seed$TACAIn+1),log(seed$AVBAIn+1), cex=sz, col=cols[cl],pch=19, xlim=c(0,5), ylim=c(-0.5,4))
image(x=log(seed$TACAIn+1),y=log(seed$AVBAIn+1),z=seeds$TACAOut )
lsds <- log(seed+1)
plot(lsds$TACAIn,lsds$TACAOut,col=mcol)
plot(lsds$TACAIn,lsds$AVBAOut,col=acol)

plot(lsds$AVBAIn,lsds$TACAOut)

#plot model estimates versus true values
#unique pairs of seedsIn and maybe mean of actual seeds out in those pairs and difference with prediction using model parameter estimates.


#differences when using a mean instead of distribution for ets seeds out bg?
#or maybe it should be treatment specific... 
#make function for all this 
#or maybe I need to make the neighbourhood smaller...













