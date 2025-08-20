plotmovement <- function(){
	par(mar=c(0.5,5,2,0))
	plot(3:timesteps,ylim=c(-0.1,1.1),type="n", xlim=c(t1+.3,timesteps), bty="n",xaxt="n", yaxt="n", ylab='', xlab='')
	x <- rep(t1,4)
	mtext("time=",1,at=t1,line=-.75,adj=1)
	title(ylab="patch",line=3.75)
	patchlab <- apply(unique(dif$origin[,c("i","j")]),1, paste,collapse=',')
	text(x=t1,y=y,paste0("(",patchlab,")"),font=2,adj=1)
	mtext(paste(parms,pnames,sep="=",collapse=" "), line=0.8, adj=0.95)
	for (t in t1:timesteps){
	#x <- rep(t,4)
		dif <- biased_diffuse(sp2[,,"seeds",t],litdifr[,,t],D=dd[t],beta=5, origins=T)
		if(dd[t]==0){x <-x+0.3} else {
			x <-x+1
			if(dd[t-1]==0){x <-x-0.7}
			out <- trackOr(dif)
			arrows(x0=x+.1,y0=out$y0, y1=out$y1-.03*out$yd/abs(out$yd), x1=x+.9, length=0.1, lwd=out$weight,col="grey")
			Pout <- ceiling(out$Pout*5)
			pcols <- colPout(Pout,ifelse(t==timesteps,T,F))
			if (percent==F){
				mapply(text, x=x[1]+exp(1/abs(out$yd)*.08)/4.25,y=out$yhlf+(out$y0-out$y1)/3, labels=round(out$out,5), cex=0.7, srt=out$srt2,col=pcols[Pout])}#out$srt)#,col='white'"red")
			if (percent==T){
				mapply(text, x=x[1]+exp(1/abs(out$yd)*.08)/4.25,y=out$yhlf+(out$y0-out$y1)/3, labels=paste(round(out$Pout*100,5),"%"), srt=out$srt2, cex=0.5,col=pcols[Pout])#col="blue"
			}
		}
		points(x,y,pch=21,cex=5, bg=hcl.colors(4, "greens", 0.5, T)[litdifr[,,t]])
		points(x,y,pch=as.character(litdifr[,,t]),cex=1.8)
		text(x,y+0.13, round(sp2[,,"inds",t],3))
		text(x,y+0.06, round(sp2[,,"seeds",t],3))
		arrows(x0=x,y0=y+0.11,y1=y+0.075,length=0.08,lwd=0.7)
		text(x+0.1,y, round(dif$N-dif$stay,3),adj=0)
		mtext(t,1,at=x,line=-.75)
	}
}
#text(x=x[f]+.5,y=((y0+y1)*.5)[f], round(na.omit(orgins$Nout),6)[f],srt=90-atan((y1-y0)[f]*0.8)/pi*180-180,cex=1.2)
	#arrows(x0=x*1.01,y0=orgins$y[1], y1=xy[xy$i==orgins[1,]$out2_i & xy$j==orgins[1,]$out2_j,]$y *.95, x1=x*1.09, length=0.1, lwd=orgins$Nout[1]*2*10^5)
#plot(1:75,pch=19,col=hcl.colors(75,"Blues 2"),cex=5)
ij2y <- function(i,j,y){
	xy <- expand.grid(i=1:nr,j=1:nc)
	xy <- matrix(1:4,nc=nc)
	z<-diag(xy[i,j])
	#xy$z <- 1:4
	#xy$y <- 	y
	out<-y[z]
	#ij <- which(xy$j%in% j & xy$i%in% i)
	#out <- xy[ij,]$y
	return(out)
}
	
trackOr <- function(dif){
	orgins <- dif$origins
	orgins$y <- rep(y,8)
	orgins$outflowT <- c(dif$N-dif$stay)
	orgins$pOut <- orgins$outflow[rep(1:4,8)]
	y0=orgins$y[!is.na(orgins$Nout)]
	outflowT <- orgins$outflowT[!is.na(orgins$Nout)]
	y1=c(na.omit(ij2y(orgins$out2_i,orgins$out2_j,y)))
	out=c(na.omit(orgins$Nout))
	out <- data.frame(y0=y0,y1=y1,out=out,outT=outflowT)
	out$Pout <- out$out/outflowT
	out$yd <- yd <- y1-y0
	out$yhlf <- (y0+y1)/2
	out$weight <- exp(out$out*5)*3-2.5#abs(log(out$out))/4 
	out$srt <- 90-atan((y1*.95-y0)*1.25)/pi*180
	out$srt2 <- yd/abs(yd)*exp(abs(yd)*.75)*35#(y1-y0)*150
	return(out)
}	

colPout <- function(Pout,colorbar=T){
	#n<-length(unique(round(Pout*100)))
	#n <- n+3
	n <- 5
	pcols <- hcl.colors(n+6,"Blues 3",rev=T)[-1:-5][-(n+(1))]
	#plot(1:n,pch=19,col=pcols,cex=5)
	#plot(1:13,pch=19,col=hcl.colors(13,"Blues 2"),cex=5)
	if (colorbar==T){
		cbx <- timesteps-((n*5):1)*.025
		segments(x0=cbx, y0=-.15,y1=-.1,col=rep(pcols,each=5),lwd=9,lend=3)
		text(cbx[c(1,(1:5)*5)],-.08,labels=c(0,(1:5)*20))
	}
	return(pcols)
}
#colPout(Pout)
#plotmovement()



