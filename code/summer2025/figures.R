#for lab retreat Aug 27 2025

#patch level Population trajectories with 
	#1)no litter no dispersal
	#2)no litter with dispersal
	#3)litter no dispersal
	#4)litter and dispersal
	
#global population trajectories x4
	
	
#plot
#pdf("../../figures/germinantDiffusion__d0.5..pdf",width=9,height=7.5)
ncol <- 6;alpha_col <- 0.9;alpha_col0<-0.5
#cols1 <- hcl.colors(ncol,"PuBu",alpha_col,rev=F)
colargs <- expand.grid(palette=c(a="PuBu",b="OrRd"),alpha=c(y=alpha_col,n=alpha_col0))
cols <- data.frame(mapply(hcl.colors, palette=colargs$palette,
               alpha=colargs$alpha,MoreArgs=list(n=ncol)))
colnames(cols) <- c(1,2,"10","20")
#cols2 <- hcl.colors(ncol,"OrRd",alpha_col,rev=F)
#col2rgb("darkgreen")
gr <- rgb(0,100/255,0,0.15)

pdf("../../figures/poppatcht-1.pdf",width=8)
par(mgp=c(2,0.1,0),mfrow=c(gridsize,gridsize),xpd=F,cex=.8,tcl=-.15,oma=c(1,.5,.5,0.1), mar=c(.8,1.3,.5,0))
t=20
for (i in 1:gridsize){
  for (j in 1:gridsize){
    plot(1:t,seq(0,1,length.out=t), type="n", ylab='',xlab='', xlim=c(1,t), ylim=c(0,0.8)) 
    lines(1:t, sp100[i,j,"inds",1:t],type='b',col=cols$"1"[1], cex=.7,lwd=2)
    lines(1:t, sp200[i,j,"inds",1:t],type='b', col=cols$"2"[1], cex=.7,lwd=2)    
  }
}
dev.off()

pdf("../../figures/poppatcht-2.pdf",width=8)
par(mgp=c(2,0.1,0),mfrow=c(gridsize,gridsize),xpd=F,cex=.8,tcl=-.15,oma=c(1,.5,.5,0.1), mar=c(.8,1.3,.5,0))
t=20
for (i in 1:gridsize){
  for (j in 1:gridsize){
    plot(1:t,seq(0,1,length.out=t), type="n", ylab='',xlab='', xlim=c(1,t), ylim=c(0,0.8)) 
    lines(1:t, sp10[i,j,"inds",1:t],type='b',col=cols$"1"[1], cex=.7,lwd=2)
    lines(1:t, sp20[i,j,"inds",1:t],type='b', col=cols$"2"[1], cex=.7,lwd=2)    
  }
}
dev.off()

pdf("../../figures/poppatcht-1and2.pdf",width=8)
par(mgp=c(2,0.1,0),mfrow=c(gridsize,gridsize),xpd=F,cex=.8,tcl=-.15,oma=c(1,.5,.5,0.1), mar=c(.8,1.3,.5,0))
t=20
for (i in 1:gridsize){
  for (j in 1:gridsize){
  	plot(1:t,seq(0,1,length.out=t), type="n", ylab='',xlab='', xlim=c(1,t), ylim=c(0,0.8)) 
    lines(1:t, sp100[i,j,"inds",1:t],type='b',col=cols$"1"[2], cex=.7,lwd=2)
    lines(1:t, sp200[i,j,"inds",1:t],type='b', col=cols$"2"[2], cex=.7,lwd=2)  
    lines(1:t, sp10[i,j,"inds",1:t],type='b',col=cols$"1"[1], cex=.7,lwd=2)
    lines(1:t, sp20[i,j,"inds",1:t],type='b',col=cols$"2"[1], cex=.7,lwd=2)    
  }
  if (i==1){
  	  legend("topright",legend=c("sp1 no dispersal","sp2 no dispersal","sp1 w dispersal","sp2 w dispersal"),bty="n",pch=1,lwd=2,col=c(cols$'1'[2],cols$'2'[2],cols$'1'[1],cols$'2'[1]),cex=.8)
  }
}
dev.off()

pdf("../../figures/poppatcht-3.pdf",width=8)
par(mgp=c(2,0.1,0),mfrow=c(gridsize,gridsize),xpd=F,cex=.8,tcl=-.15,oma=c(1,.5,.5,0.1), mar=c(.8,1.3,.5,0))
t=20
for (i in 1:gridsize){
  for (j in 1:gridsize){
  	plot(1:t,seq(0,1,length.out=t), type="n", ylab='',xlab='', xlim=c(1,t), ylim=c(0,0.8))  
  	lines(c(0,rep(1:t,each=2))+0.5,c(rep(lit[i,j,1:t], each=2),lit[i,j,t+1]),col="lightgrey",lwd=2,type='l')
    lines(1:t, sp10l[i,j,"inds",1:t],type='b',col=cols$"1"[1], cex=.7,lwd=2)
    lines(1:t, sp20l[i,j,"inds",1:t],type='b',col=cols$"2"[1], cex=.7,lwd=2)    
  }
}
dev.off()

pdf("../../figures/poppatcht-4.pdf",width=8)
par(mgp=c(2,0.1,0),mfrow=c(gridsize,gridsize),xpd=F,cex=.8,tcl=-.15,oma=c(1,.5,.5,0.1), mar=c(.8,1.3,.5,0))
t=20
for (i in 1:gridsize){
  for (j in 1:gridsize){
  	plot(1:t,seq(0,1,length.out=t), type="n", ylab='',xlab='', xlim=c(1,t), ylim=c(0,0.8))  
  	lines(c(0,rep(1:t,each=2))+0.5,c(rep(lit[i,j,1:t], each=2),lit[i,j,t+1]),col="lightgrey",lwd=2,type='l')
    lines(1:t, sp1[i,j,"inds",1:t],type='b',col=cols$"1"[1], cex=.7,lwd=2)
    lines(1:t, sp2[i,j,"inds",1:t],type='b',col=cols$"2"[1], cex=.7,lwd=2)    
  }
}
dev.off()

pdf("../../figures/poppatcht-234.pdf",width=8)
par(mgp=c(2,0.1,0),mfrow=c(gridsize,gridsize),xpd=F,cex=.8,tcl=-.15,oma=c(1,.5,.5,0.1), mar=c(.8,1.3,.5,0))
t=20
for (i in 1:gridsize){
  for (j in 1:gridsize){
  	plot(1:t,seq(0,1,length.out=t), type="n", ylab='',xlab='', xlim=c(1,t), ylim=c(0,0.8))  
  	lines(c(0,rep(1:t,each=2))+0.5,c(rep(lit[i,j,1:t], each=2),lit[i,j,t+1]),col="lightgrey",lwd=2,type='l')
  	
  	lines(1:t, sp10[i,j,"inds",1:t],type='b',col=cols$"1"[3], cex=.7,lwd=2)
    lines(1:t, sp20[i,j,"inds",1:t],type='b',col=cols$"2"[3], cex=.7,lwd=2) 
  	lines(1:t, sp10l[i,j,"inds",1:t],type='b',col=cols$"1"[2], cex=.7,lwd=2)
    lines(1:t, sp20l[i,j,"inds",1:t],type='b',col=cols$"2"[2], cex=.7,lwd=2)
    lines(1:t, sp1[i,j,"inds",1:t],type='b',col=cols$"1"[1], cex=.7,lwd=2)
    lines(1:t, sp2[i,j,"inds",1:t],type='b',col=cols$"2"[1], cex=.7,lwd=2)    
  }
  if (i==1){
  	legend("topright",legend=paste(rep(c("sp1","sp2"),each=3),rep(c("dispersal no lit","lit no dispersal", "litter and dispersal"),2)) ,bty="n", pch=1,lwd=2, col=c(cols$'1'[3:1],cols$'2'[3:1]),cex=.8)
  }
}
dev.off()

pdf("../../figures/popglob.pdf",width=11)
par(mgp=c(2,0.1,0),mfrow=c(1,2),xpd=F,cex=.8,tcl=-.15,oma=c(1,.5,.5,0.1), mar=c(.8,1.3,.5,0))
t=20
fun <- mean
plot(1:t,seq(0,1,length.out=t), type="n", ylab='',xlab='', xlim=c(1,t), ylim=c(0,0.5))  
litsum <- apply(lit,3,fun)
lines(c(0,rep(1:t,each=2))+0.5,c(rep(litsum[1:t], each=2),litsum[t+1]),col="lightgrey",lwd=2,type='l')
lines(1:t, apply(sp1[,,"inds",], 3,fun)[1:t], type='b', col=cols$"1"[1], cex=.7,lwd=2)
lines(1:t, apply(sp2[,,"inds",], 3,fun)[1:t], type='b', col=cols$"2"[1], cex=.7,lwd=2)   
lines(1:t, apply(sp10[,,"inds",], 3,fun)[1:t], type='b', col=cols$"1"[3], cex=.7,lwd=2)
lines(1:t, apply(sp20[,,"inds",], 3,fun)[1:t], type='b', col=cols$"2"[3], cex=.7,lwd=2)    
lines(1:t, apply(sp10l[,,"inds",], 3,fun)[1:t], type='b', col=cols$"1"[2], cex=.7,lwd=2)
lines(1:t, apply(sp20l[,,"inds",], 3,fun)[1:t], type='b', col=cols$"2"[2], cex=.7,lwd=2)  
legend("bottomleft",legend=c("dispersal no lit","lit no dispersal","litter and dispersal"), bty="n",pch=1, lwd=2, col=c(cols$'1'[3],cols$'1'[2],cols$'1'[1]),cex=1, title="sp1")
legend("topleft",legend=c("dispersal no lit","lit no dispersal","litter and dispersal"), bty="n",pch=1, lwd=2, col=c(cols$'2'[3:1]),cex=1, title="sp2")
fun <- sum
plot(1:t,seq(0,1,length.out=t), type="n", ylab='',xlab='', xlim=c(1,t), ylim=c(0,2))  
litsum <- apply(lit,3,fun)
lines(c(0,rep(1:t,each=2))+0.5,c(rep(litsum[1:t], each=2),litsum[t+1]),col="lightgrey",lwd=2,type='l')
lines(1:t, apply(sp1[,,"inds",], 3,fun)[1:t], type='b', col=cols$"1"[1], cex=.7,lwd=2)
lines(1:t, apply(sp2[,,"inds",], 3,fun)[1:t], type='b', col=cols$"2"[1], cex=.7,lwd=2)   
lines(1:t, apply(sp10[,,"inds",], 3,fun)[1:t], type='b', col=cols$"1"[3], cex=.7,lwd=2)
lines(1:t, apply(sp20[,,"inds",], 3,fun)[1:t], type='b', col=cols$"2"[3], cex=.7,lwd=2)    
lines(1:t, apply(sp10l[,,"inds",], 3,fun)[1:t], type='b', col=cols$"1"[2], cex=.7,lwd=2)
lines(1:t, apply(sp20l[,,"inds",], 3,fun)[1:t], type='b', col=cols$"2"[2], cex=.7,lwd=2)  
 
dev.off()

saveRDS(data.frame(parms,pnames),"parms_figs1.RDS")





