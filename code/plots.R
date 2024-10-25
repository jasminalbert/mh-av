datLoc <- "../data/"
filename <- "MHAVseeds.RDS"
mas <- readRDS(paste0(datLoc,filename))
seed<-readRDS("si-so.RDS")
modelOut <-readRDS("modelest1025.RDS")



table(seed$TACAIn)
tab <- table(seed$AVBAIn,seed$TACAIn)
ain <- sort(unique(seed$AVBAIn))
min <- sort(unique(seed$TACAIn))
bgc <- rgb(0,0,0,0.02)
ptc <- rgb(0,0,0,0.1)
ins <- cbind(ain,min)
nuobs <- sum(tab>0)
uobs <- which(tab>0,arr.ind=T)
plot(seq_along(ain),seq_along(min),xaxt="n",yaxt="n",xlab="av",ylab="mh",type="n");axis(1,labels=round(ain),at=seq_along(ain));axis(2,labels=round(min),at=seq_along(min))

for (n in 1:nuobs){
	atmp <- round(as.numeric(rownames(tab)[uobs[n,"row"]]),1)
	mtmp <- round(as.numeric(colnames(tab)[uobs[n,"col"]]),1)
	print(cbind(atmp,mtmp))
	seedtmp <- seed[round(seed$TACAIn,1)==mtmp& round(seed$AVBAIn,1)==atmp,]
	tou <- seedtmp$TACAOut;aou<-seedtmp$AVBAOut
	tost <- c(mn=mean(tou,na.rm=T),sd=sd(tou,na.rm=T))
	aost <- c(mn=mean(aou,na.rm=T),sd=sd(aou,na.rm=T))
	print(aost)
	#hist(tou);hist(aou)
	psa <- round(log(aou+1))+1
	psm <- round(log(tou+1))+1
	x <- which(round(ain,1)==atmp);y <- which(round(min,1)==mtmp)
	x <- x[rep(1,length(aou))]; y <- y[rep(1,length(aou))]
	points(x,y,cex=psa,pch=19,col=bgc,bg=bgc)
}
plot(seq_along(ain),seq_along(min),xaxt="n",yaxt="n",xlab="av",ylab="mh",type="n");axis(1,labels=round(ain),at=seq_along(ain));axis(2,labels=round(min),at=seq_along(min))

for (n in 1:nuobs){
	atmp <- round(as.numeric(rownames(tab)[uobs[n,"row"]]),1)
	mtmp <- round(as.numeric(colnames(tab)[uobs[n,"col"]]),1)
	print(cbind(atmp,mtmp))
	seedtmp <- seed[round(seed$TACAIn,1)==mtmp& round(seed$AVBAIn,1)==atmp,]
	tou <- seedtmp$TACAOut;aou<-seedtmp$AVBAOut
	tost <- c(mn=mean(tou,na.rm=T),sd=sd(tou,na.rm=T))
	aost <- c(mn=mean(aou,na.rm=T),sd=sd(aou,na.rm=T))
	print(aost)
	#hist(tou);hist(aou)
	psa <- round(log(aou+1))+1
	psm <- round(log(tou+1))+1
	x <- which(round(ain,1)==atmp);y <- which(round(min,1)==mtmp)
	x <- x[rep(1,length(tou))]; y <- y[rep(1,length(tou))]
	points(x,y,cex=psm,pch=19,col=bgc,bg=bgc)
}
#real value plots now simulate predicted values and check difference
#instead f 4x4 could just have x axis be 16 differnt categories






