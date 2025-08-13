#working on (1) burn in monoculture
#storage array
dimnames <- list(i=1:gridsize, j=1:gridsize,stage=c("inds","seeds","dispered"), 
                 time=1:(timesteps+1))
sp1 <- array(NA, dim=c(gridsize,gridsize,3,timesteps+1), dimnames)
sp2 <- sp1
sp1[,,"inds",1] <- N1
sp2[,,"inds",1] <- 0
sp2[,,"inds",2] <- N2
#litter
lit <- sp1[,,1,]
lit[,,1] <- ifelse(is.na(N1),NA,0)

#no litter models
sp10 <- sp1; sp20 <- sp2
t=1
for (t in 1:timesteps){
  #seed=sample(1:1e6,1)
  litdif <- litDif(lit[,,t])
  #litdifnorm <- litdif/sum(abs(litdif)) 
  lambda1lit <- lambda1*exp(-litdif*rho)#no lit at t=1
  sp1[,,"seeds",t] <- lambda1lit*sp1[,,"inds",t]*(1-sp1[,,"inds",t]-alpha*sp2[,,"inds",t])
  sp10[,,"seeds",t] <- lambda1*sp10[,,"inds",t]*(1-sp10[,,"inds",t]-alpha*sp20[,,"inds",t])
  sp2[,,"seeds",t] <- lambda2*sp2[,,"inds",t]*(1-sp2[,,"inds",t]-beta*sp1[,,"inds",t])
  sp20[,,"seeds",t] <- lambda2*sp20[,,"inds",t]*(1-sp20[,,"inds",t]-beta*sp10[,,"inds",t])
  zero <- mapply(function(x){x[,,"seeds",t]},list(sp1,sp10,sp2,sp20))<0
  if (any(zero)){cat("\nZERO!",t)}
  #mapply(function(x){x[,,"seeds",t]},list(sp1,sp10,sp2,sp20))
  if(any(sp1[,,"seeds",t] <0)){sp1[,,"seeds",t][sp1[,,"seeds",t]<0]<-0}
  if(any(sp2[,,"seeds",t] <0)){sp2[,,"seeds",t][sp2[,,"seeds",t]<0]<-0}
  #dispersal

  litdif[,] <- rank(-litdif)
  sp1[,,"dispered",t] <- biased_diffuse(sp1[,,"seeds",t],lit[,,t],D=D,beta=0)$N_new
  sp10[,,"dispered",t] <- biased_diffuse(sp10[,,"seeds",t],lit[,,t],D=D,beta=0)$N_new
  sp2[,,"dispered",t] <- biased_diffuse(sp2[,,"seeds",t],litdif,D=D,beta=5)$N_new
  sp20[,,"dispered",t] <- biased_diffuse(sp20[,,"seeds",t],lit[,,t],D=D,beta=0)$N_new
  #sp2[,,"dispered",t] <- difdis(D*theta*litdif, sp2[,,"seeds",t])
  #deposit litter
  lit1 <- sp1[,,"inds",t]*l1
  lit2 <- sp2[,,"inds",t]*l2
  lit[,,t+1] <- lit1+lit2
  #t+1 inds are what dispersed
  sp1[,,"inds",t+1] <- sp1[,,"dispered",t]
  sp10[,,"inds",t+1] <- sp10[,,"dispered",t]
  if (t>1){
    sp2[,,"inds",t+1] <- sp2[,,"dispered",t]
    sp20[,,"inds",t+1] <- sp20[,,"dispered",t]
  }
}
litmaxes <- apply(lit[,,-1], 3, function(x) which(x==max(x),T ))
litmax <- lit
litmax[,,] <- FALSE
litmax[cbind(t(litmaxes),colnames(litmaxes))] <- TRUE

##mapply(function(x){x<0}, list(sp1[,,"dispered",],sp10[,,"dispered",],sp2[,,"dispered",],sp20[,,"dispered",] ))
#plot
#pdf("../../figures/germinantDiffusion__d0.5..pdf",width=9,height=7.5)
pnames <- c(lambda1,lambda2,alpha,beta,rho,theta,l1,l2,D)
par(mgp=c(2,0.1,0),mfrow=c(2,2),xpd=F,cex=.8,tcl=-.15,oma=c(1,.5,.5,0), mar=c(.8,1.3,.5,0))
for (i in 1:2){
  for (j in 1:2){
    plot(1:t,seq(0,1,length.out=t),type="n",ylab='',xlab='',xlim=c(1,t),
         ylim=c(0,0.8))
    polygon(x=c(1:t-.3,t*1.02,t*1.02,t:1-.3),y=c(rep(-0.05,(t+1)),rev(lit[i,j,1:t])[1],rev(lit[i,j,1:t])),
            col=gr,border=NA)
    lines(1:t-.3,lit[i,j,1:t],col=ifelse(litmax[i,j,1:t],"limegreen","darkgreen"),type='b',pch=ifelse(litmax[i,j,1:t],8,3),
          cex=ifelse(litmax[i,j,1:t],1.5,0.6),lwd=0.4)
    lines(1:t+.5,sp1[i,j,"dispered",1:t],type='b',col=cols1[3],lwd=0.5,pch=17)
    lines(2:t+.5,sp2[i,j,"dispered",2:t],type='b',col=cols2[3],lwd=0.5,pch=17)
    lines(1:t+.5,sp10[i,j,"dispered",1:t],type='b',col=cols1[3],lwd=0.5,pch=2,lty=3,cex=0.7)
    lines(2:t+.5,sp20[i,j,"dispered",2:t],type='b',col=cols2[3],lwd=0.5,pch=2,lty=3,cex=0.7)
    lines(sp1[i,j,"inds",1:t],type='b',col=cols1[1])
    lines(2:t,sp2[i,j,"inds",2:t],type='b',col=cols2[1])
    lines(sp10[i,j,"inds",1:t],type='b',lty=2,col=cols1[1])
    lines(2:t,sp20[i,j,"inds",2:t],type='b',col=cols2[1],lty=2)
    points(1:t+.2,sp1[i,j,"seeds",1:t],pch=18,col=cols1[2])
    points(2:t+.2,sp2[i,j,"seeds",2:t],pch=18,col=cols2[2])
    points(1:t+.2,sp10[i,j,"seeds",1:t],pch=5,col=cols1[2],cex=0.7)
    points(2:t+.2,sp20[i,j,"seeds",2:t],pch=5,col=cols2[2],cex=0.7)
    
    #legend("topleft", legend=c("sp1","sp1notlit","seedslit","seedsnolit",
    #                               "disp_lit","dispNolit","lit"),col=c(1,"darkred","blue","orange",
    #                                                                  "skyblue","tan",3),lty=1,bty="n")
  }
}
title(main=paste(parms,pnames,sep="=",collapse=" "),outer=T,font.main=1,line=-.3)
#dev.off()
