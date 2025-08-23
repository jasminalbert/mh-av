plotmovement <- function(sparray,theta,...){
	par(mar=c(0.5,5,2,0))
	plot(t1:timesteps,ylim=c(-0.1,1.1),type="n", xlim=c(t1+.3,timesteps), bty="n",xaxt="n", yaxt="n", ylab='', xlab='')
	x <- rep(t1,gridsize^2)
	mtext("time=",1,at=t1,line=-.75,adj=1)
	title(ylab="patch",line=3.75)
	mtext(paste(parms,pnames,sep="=",collapse=" "), line=0.8, adj=0.95)
	for (t in t1:timesteps){
	#x <- rep(t,4)
		dif <- biased_diffuse(sparray[,,"seeds",t],litdifr[,,t],D=dd[t],beta=theta, origins=T)
		cat(t)
		if(dd[t]==0){x <-x+0.3} else {
			x <-x+1
			if(dd[t-1]==0){x <-x-0.7}
			out <- trackOr(dif)
			arrows(x0=x+.1,y0=out$y0, y1=out$y1-.03*out$yd/abs(out$yd), x1=x+.9, length=0.1, lwd=out$weight,col="grey")
			Pout <- ceiling(out$Pout*100)
			pcols <- colPout(Pout,ifelse(t==timesteps,T,F),...)
			if (percent==F){
				mapply(text, x=x[1]+exp(1/abs(out$yd)*.08)/4.25,y=out$yhlf+(out$y0-out$y1)/3, labels=round(out$out,5), cex=0.8, srt=out$srt2,col=pcols[Pout])}#out$srt)#,col='white'"red")
			if (percent==T){
				mapply(text, x=x[1]+exp(1/abs(out$yd)*.08)/4.25,y=out$yhlf+(out$y0-out$y1)/3, labels=paste(round(out$Pout*100,5),"%"), srt=out$srt2, cex=0.8,col=pcols[Pout])#col="blue"
			}
		}
		points(x,y,pch=21,cex=5, bg=hcl.colors(gridsize^2, "greens", 0.5, T)[litdifr[,,t-1]])
		points(x,y,pch=as.character(litdifr[,,t-1]),cex=1.8)
		text(x,y+0.13, round(sparray[,,"inds",t],3))
		text(x,y+0.06, round(sparray[,,"seeds",t],3))
		pal <- colorRampPalette(c("red","grey","green"))(15)
		scal <- seq(0,1,length.out=15)
		seedssign <- data.frame(scal,pal)
		pal[1 + floor(scal*(15-1))]
		arrows(x0=x,y0=y+0.11,y1=y+0.076,length=0.06,lwd=1.2,lend=2, col= seedssign[round(c(sparray[,,"seeds",t]/sparray[,,"inds",t])*10),2] )
		text(x+0.1,y, round(dif$N-dif$stay,3),adj=0)
		mtext(t,1,at=x,line=-.75)
	}
	patchlab <- apply(dif$origin[,c("i","j"),1],1, paste,collapse=',')
	text(x=t1,y=y,paste0("(",patchlab,")"),font=2,adj=1)
}
#text(x=x[f]+.5,y=((y0+y1)*.5)[f], round(na.omit(orgins$Nout),6)[f],srt=90-atan((y1-y0)[f]*0.8)/pi*180-180,cex=1.2)
	#arrows(x0=x*1.01,y0=orgins$y[1], y1=xy[xy$i==orgins[1,]$out2_i & xy$j==orgins[1,]$out2_j,]$y *.95, x1=x*1.09, length=0.1, lwd=orgins$Nout[1]*2*10^5)
#plot(1:75,pch=19,col=hcl.colors(75,"Blues 2"),cex=5)
ij2y <- function(i,j,y){
	xy <- expand.grid(i=1:gridsize,j=1:gridsize)
	xy <- matrix(1:gridsize^2,nc=gridsize)
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
	ynbs <- rep(y,n_nbs)
	outflowT <- c(dif$N-dif$stay)
	pOut <- outflowT[rep(1:gridsize^2,n_nbs)]
	y0=ynbs[!is.na(orgins[,"Nout",])]
	outflowT <- pOut[!is.na(orgins[,"Nout",])]
	y1=c(na.omit(ij2y(orgins[,"out2_i",],orgins[,"out2_j",],y)))
	out=c(na.omit(c(orgins[,"Nout",])))
	out <- data.frame(y0=y0,y1=y1,out=out,outT=outflowT)
	out$Pout <- out$out/outflowT
	out$yd <- yd <- y1-y0
	out$yhlf <- (y0+y1)/2
	out$weight <- exp(out$out*15)*0.5#abs(log(out$out))/4 
	#*5)*3-2.5
	out$srt <- 90-atan((y1*.95-y0)*1.25)/pi*180
	out$srt2 <- yd/abs(yd)*exp(abs(yd)*.75)*35#(y1-y0)*150
	return(out)
}	
stretch_high <- function(x, power = 2) {
  ( (x - min(x)) / (max(x) - min(x)) )^power
}
colPout <- function(Pout,colorbar=T,pall="Viridis", trim=0,w){
	#n<-length(unique(round(Pout*100)))
	#n <- n+3
	n<-5
	vals <- 1:100
	ncol <- length(vals)
	pal <- hcl.colors(ncol,pall,rev=T)
	if (trim>0){
		pcols <- hcl.colors(n+trim,pall,rev=T)[-1:-trim]
	}
	scaled <- stretch_high(vals, power = w)
	pcols   <- pal[1 + floor(scaled*(ncol-1))]
	#plot(1:n,pch=19,col=pcols,cex=5)
	#plot(1:13,pch=19,col=hcl.colors(13,"Blues 2"),cex=5)
	if (colorbar==T){
		cbx <- t1+(seq(-1,1.5,length.out=ncol))*.5
		segments(x0=cbx, y0=1.18,y1=1.13,col=pcols,lwd=9,lend=3)
		text(cbx[c(1,(1:5)*20)],1.2,labels=c(0,(1:5)*20))
	}
	return(pcols)
}
#axis(1)
#colPout(Pout)
#plotmovement()


#plot(vals,stretch_high(vals,3),col=rainbow(ncol,0.7,0.9),pch=19,cex=5)
#points(vals,log(vals)/max(log(vals)),col=rainbow(ncol,0.7,0.9),pch=19,cex=5)
