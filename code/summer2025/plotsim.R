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

pnames <- c(lambda1,lambda2,alpha,beta,rho,theta,l1,l2,D,LLfun)
par(mgp=c(2,0.1,0),mfrow=c(2,2),xpd=F,cex=.8,tcl=-.15,oma=c(1,.5,.5,0), mar=c(.8,1.3,.5,0))
for (i in 1:2){
  for (j in 1:2){
    plot(1:t,seq(0,1,length.out=t),type="n",ylab='',xlab='',xlim=c(1,t),
         ylim=c(0,0.8))
    polygon(x=c(1:t-.3,t*1.02,t*1.02,t:1-.3),y=c(rep(-0.05,(t+1)),rev(lit[i,j,1:t])[1],rev(lit[i,j,1:t])),
            col=gr,border=NA)
    lines(1:t-.3,lit[i,j,1:t],col=ifelse(litmax[i,j,1:t],"limegreen","darkgreen"),type='b',pch=ifelse(litmax[i,j,1:t],8,3),
          cex=ifelse(litmax[i,j,1:t],1.5,0.6),lwd=0.4)
    if(plot_equil==T){
    	lambda1 <- lambda[i,j,-(t+1)]
    	points(equils()$xstary0,pch="--",col=cols$'10'[2],cex=2)
    	points(rep(equils()$ystarx0,t),pch="--",col=cols$'20'[2],cex=2)
    	points(equils()$xstar,pch="--",col=cols$'1'[1],cex=2)
    	points(equils()$ystar,pch="--",col=cols$'2'[1],cex=2)
    }
    if(plot_zngi==T){
    	lambda1 <- lambda[i,j,-(t+1)]
    	points(isocline_y(sp1[i,j,"inds",-(t+1)]),pch="-",col="magenta",cex=2)
    	points(isocline_x(sp2[i,j,"inds",-(t+1)]),pch="-",col="cyan",cex=2)
    }
    if(plot_dispersed==T){
    lines(1:t+.5,sp1[i,j,"dispered",1:t],type='b',col=cols$"1"[3],lwd=0.5,pch=17)
    lines(1:t+.5,sp2[i,j,"dispered",1:t],type='b',col=cols$"2"[3],lwd=0.5,pch=17)
    lines(1:t+.5,sp10[i,j,"dispered",1:t],type='b',col=cols$"10"[3],lwd=0.5,pch=2,lty=3,cex=0.7)
    lines(1:t+.5,sp20[i,j,"dispered",1:t],type='b',col=cols$"20"[3],lwd=0.5,pch=2,lty=3,cex=0.7)}
    
    lines(sp1[i,j,"inds",1:t],type='b',col=cols$"1"[1],cex=.7)
    lines(1:t,sp2[i,j,"inds",1:t],type='b',col=cols$"2"[1],cex=.7)
    lines(sp10[i,j,"inds",1:t],type='b',lty=2,col=cols$"10"[1],cex=.7)
    lines(1:t,sp20[i,j,"inds",1:t],type='b',col=cols$"20"[1],lty=2,cex=.7)
    points(1:t+.2,sp1[i,j,"seeds",1:t],pch=18,col=cols$"1"[2])
    points(1:t+.2,sp2[i,j,"seeds",1:t],pch=18,col=cols$"2"[2])
    points(1:t+.2,sp10[i,j,"seeds",1:t],pch=5,col=cols$"10"[2],cex=0.7)
    points(1:t+.2,sp20[i,j,"seeds",1:t],pch=5,col=cols$"20"[2],cex=0.7)
    
    #legend("topleft", legend=c("sp1","sp1notlit","seedslit","seedsnolit",
    #                               "disp_lit","dispNolit","lit"),col=c(1,"darkred","blue","orange",
    #                                                                  "skyblue","tan",3),lty=1,bty="n")
  }
}
title(main=paste(parms,pnames,sep="=",collapse=" "),outer=T,font.main=1,line=-.3)
#dev.off()
