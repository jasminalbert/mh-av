
#From meeting on 7 23 25 w LH:
#  Talked about using diffusion on a germinant model and movement of germinants is determined by discrete diffusion
#Using a model with reduced parameters by scaling of lambda and alphas
#Dont accumulate litter just make litter a fraction of previous years counts

#scaled parameters

#lambda = lambda2/lambda1
#alpha = alpha12/alpha11
#beta = alpha21/alpha22

#without dispersal
#N_1(t+1) = N_1(t)(1-N_1(t)-alpha*N_2(t))
#N_2(t+1) = lambda*N_2(t)(1-N_2(t)-beta*N_1(t))

#population growth function
growth <- function(lambda1, lambda2, alpha, beta, N){
  N1 <- N[1]
  N2 <- N[2]
  N1next = lambda1*N1*(1-N1-alpha*N2)
  N2next = lambda2*N2*(1-N2-beta*N1)
  return(c(N1next,N2next))
}
lambda1 <- lambda2 <- 1
alpha <- 0.1
beta <- 0.2
growth(lambda=1,alpha,beta,N)
(1-1/1.5)/beta
alpha <- beta12; beta <- beta21
timesteps <- 20
N <- array(c(N1,N2),dim=c(dim(N1),2))
N1 = matrix(0.6)
N2 = matrix(0.5)
N <- array(NA, dim=c(dim(N1)[1],timesteps,2))
N[,1,] <- c(N1,N2) 
for (t in 2:timesteps){
  N[,t,] <- growth(lambda1,lambda2,alpha,beta,N[,t-1,])
}
par(mfrow=c(1,2),mar=c(2,1,1,1),oma=c(3,3,0,0))
plot(1:timesteps,N[,,1],type='l', ylim=range(N));lines(1:timesteps,N[,,2],col=2)
abline(h=eq['ystar'],col=2,lty=2,lwd=0.5);abline(h=eq['xstar'],lty=2,lwd=0.5)
plot(N[,,1],N[,,2],type='l',col='orange',
     xlim=range(N[,,1])*c(.7,1.3),ylim=range(N[,,2])*c(.7,1.3))
points(x=eq['xstar'],y=eq['ystar'],xpd=NA,pch=8,col=4);points(
  N[,t,1],N[,t,2],pch=19,cex=0.5);points(N[,1,1],N[,1,2],cex=0.5)


#isoclines
N1zngi <- function(alpha, N2){
  return(-alpha*N2)
}
N2zngi <- function(lambda, beta, N1){
  return(((lambda-1)/lambda)-beta*N1)
}

n <- seq(-1,1,0.01)
plot(N1zngi(alpha,n),n,type='l',xlab="N1",ylab="N2",xlim=range(n), ylim=range(n))
lines(n,N2zngi(lambda,beta,n),col=2)
abline(h=0,lty=2,col="grey")
abline(v=0,lty=2,col="grey")



lambda1 <- 1.7
lambda2 <- 1.5
alpha <- 0.4
beta <- 0.55
D <- 0.3 #diffusion coefficient
N1 <- matrix(c(0.8,0.1),ncol=2)
N2 <- matrix(c(0.1,0.8),ncol=2)
sp1array <- array(NA,dim=c(dim(N1),timesteps))
sp2array <- sp1array
sp1array[,,1] <- N1; sp2array[,,1] <- N2 #germinants at time 1
for (t in 2:timesteps){
  #germinants compete to make seeds
  N1seeds <- lambda1*N1*(1-N1-alpha*N2)
  N2seeds <- lambda2*N2*(1-N2-beta*N1)
  #litter
  N1litter <- N1
  N2litter <- 0.3 * N2
  litter <- N1litter+N2litter
  #how these seeds turn into germinants is determined by diffusion
  diff1 <- D*(sum(N1seeds)-N1seeds-N1seeds) #works for 1x2 but will have to change
  diff2 <- D*(sum(N2seeds)-N2seeds-N2seeds)
  diff1 <- D*(sum(litter)-litter-litter) #works for 1x2 but will have to change
  diff2 <- D*(sum(litter)-litter-litter)
  #sp1 germinant diffusion does not depend on litter
  #do we let it be random then?
  N1 <- N1seeds + diff1
  N2 <- N2seeds + diff2
  
  sp1array[,,t] <- N1
  sp2array[,,t] <- N2
}

for (i in 1){#edit for generalizing gridsize
  for (j in 1:2){
    plot(1:timesteps,sp1array[i,j,],type='l',ylim=c(0,1))
    lines(1:timesteps,sp2array[i,j,],col=2)
  }
}



lambda1 <- 1.7
lambda2 <- 1.5
alpha <- 0.4
beta <- 0.55
D <- 0.3 #diffusion coefficient
N1 <- matrix(c(0.8,0.1),ncol=2)
N2 <- matrix(c(0.1,0.8),ncol=2)
sp1array <- array(NA,dim=c(dim(N1),timesteps))
sp2array <- sp1array
sp1array[,,1] <- N1; sp2array[,,1] <- N2 #germinants at time 1
for (t in 2:timesteps){
  #germinants compete to make seeds
  N1seeds <- lambda1*N1*(1-N1-alpha*N2)
  N2seeds <- lambda2*N2*(1-N2-beta*N1)
  #litter
  N1litter <- N1
  N2litter <- 0.3 * N2
  litter <- N1litter+N2litter
  #how these seeds turn into germinants is determined by diffusion
  diff1 <- D*(sum(N1seeds)-N1seeds-N1seeds) #works for 1x2 but will have to change
  diff2 <- D*(sum(N2seeds)-N2seeds-N2seeds)
  diff1 <- D*(sum(litter)-litter-litter) #works for 1x2 but will have to change
  diff2 <- D*(sum(litter)-litter-litter)
  #sp1 germinant diffusion does not depend on litter
  #do we let it be random then?
  N1 <- N1seeds + diff1
  N2 <- N2seeds + diff2
  
  sp1array[,,t] <- N1
  sp2array[,,t] <- N2
}

for (i in 1){#edit for generalizing gridsize
  for (j in 1:2){
    plot(1:timesteps,sp1array[i,j,],type='l',ylim=c(0,1))
    lines(1:timesteps,sp2array[i,j,],col=2)
  }
}

####
timesteps=5
gridsize <- 2
#random dispersal array and function
rd_ar <- array(NA,dim=c(gridsize,gridsize,2,timesteps))
dimnames(rd_ar)<-list(i=1:2,j=1:2,sp=c('sp1','sp2'),t=1:timesteps)
for(s in dimnames(rd_ar)[[3]]){
  for (t in 1:timesteps){
    rd_ar[,,s,t] <- rnorm(gridsize^2)
  }
}
difdis <- function(rd,seeds){
  dispin <- rd*seeds
  dispout=(-rd*seeds)/3
  dis=array(NA,dim=dim(seeds))
  for (i in 1:(ncol(seeds)^2)){
    dis[i] = seeds[i]-sum(dispout[-i]) #whats left after you gave
    #dis_1[i] = seeds[i]-rd[i]*seeds[i]
  }
  res <- dis - dispin
  res[is.na(res)] <-0
  res1 <- res
  if (any(res<0)){
    neg <- res[res<0]
    add <- sum(abs(neg))+0.001*length(neg)
    res[!res<0] <- res[!res<0] - add/length(res[!res<0])
    res[res<0] <- 0.001 
  }
  if(abs(sum(res1)-sum(res))>1e-6){stop("sum(res1)-sum(res))>1e-6")}
  if(abs(sum(res1)-sum(seeds))>1e-6){stop("sum(res1)-sum(seeds))>1e-6")}
  return(res)
}
l1 <- 0.5;l2<-0.1
rho <- 0.5
theta <- 1
D <- 0.8
lambda1 <- 1.7
lambda2 <- 1.5
alpha <- 0.4
beta <- 0.55
#make a loop
#individuals at time 1
N1 <- matrix(c(0.8,0.05,0.1,0.7),ncol=gridsize) 
N2 <- 0.85-N1
#storage array
dimnames <- list(i=1:gridsize, j=1:gridsize,stage=c("inds","seeds","dispered"), 
                 time=1:timesteps)
sp1 <- array(NA, dim=c(gridsize,gridsize,3,timesteps), dimnames)
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
  litdifnorm <- litdif/sum(abs(litdif)) 
  lambda1lit <- lambda1*(1-rho*litdif) #no lit at t=1
  sp1[,,"seeds",t] <- lambda1lit*sp1[,,"inds",t]*(1-sp1[,,"inds",t]-alpha*sp2[,,"inds",t])
  sp10[,,"seeds",t] <- lambda1*sp10[,,"inds",t]*(1-sp10[,,"inds",t]-alpha*sp20[,,"inds",t])
  sp2[,,"seeds",t] <- lambda2*sp2[,,"inds",t]*(1-sp2[,,"inds",t]-beta*sp1[,,"inds",t])
  sp20[,,"seeds",t] <- lambda2*sp20[,,"inds",t]*(1-sp20[,,"inds",t]-beta*sp10[,,"inds",t])
  #dispersal
  randis1 <- litDif(rd_ar[,,"sp1",t])
  randis1 <- randis1/sum(abs(randis1))#sum(randis1)
  randis2 <- litDif(rd_ar[,,"sp2",t])
  randis2 <- randis2/sum(abs(randis2))#sum(randis2)
  sp1[,,"dispered",t] <- difdis(randis1*D,sp1[,,"seeds",t])
  sp10[,,"dispered",t] <- difdis(randis1*D,sp10[,,"seeds",t])
  sp2[,,"dispered",t] <- difdis(D*theta*litdifnorm, sp2[,,"seeds",t])
  sp20[,,"dispered",t] <- difdis(randis2*D,sp20[,,"seeds",t])
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
mapply(function(x){x<0}, list(sp1[,,"dispered",],sp10[,,"dispered",],sp2[,,"dispered",],sp20[,,"dispered",] ))

parms <- c("lambda1","lambda2","alpha","beta","rho","theta","l1","l2","D")
pnames <- c(lambda1,lambda2,alpha,beta,rho,theta,l1,l2,D)
ncol <- 6;alpha_col <- 0.9
cols1 <- hcl.colors(ncol,"PuBu",alpha_col,rev=F)
cols2 <- hcl.colors(ncol,"OrRd",alpha_col,rev=F)
col2rgb("darkgreen")
gr <- rgb(0,100/255,0,0.15)
i=1;j=1
theta=1
t=5

pdf("../../figures/germinantDiffusion__d0.8.pdf",width=9,height=7.5)
par(mgp=c(2,0.1,0),mfrow=c(2,2),xpd=F,cex=.8,tcl=-.15,oma=c(1,.5,.5,0), mar=c(.8,1.3,.5,0))
for (i in 1:2){
  for (j in 1:2){
    plot(1:t,seq(0,1,length.out=t),type="n",ylab='',xlab='',xlim=c(1,t),
         ylim=c(0,0.8))
    polygon(x=c(1:t-.3,t*1.02,t*1.02,t:1-.3),y=c(rep(-0.05,6),rev(lit[i,j,])[1],rev(lit[i,j,])),
            col=gr,border=NA)
    lines(1:t-.3,lit[i,j,],col="darkgreen",type='b',pch=3,cex=0.6,lwd=0.4)
    lines(1:t+.5,sp1[i,j,"dispered",],type='b',col=cols1[3],lwd=0.5,pch=17)
    lines(1:t+.5,sp2[i,j,"dispered",],type='b',col=cols2[3],lwd=0.5,pch=17)
    lines(1:t+.5,sp10[i,j,"dispered",],type='b',col=cols1[3],lwd=0.5,pch=2,lty=3,cex=0.7)
    lines(2:t+.5,sp20[i,j,"dispered",2:t],type='b',col=cols2[3],lwd=0.5,pch=2,lty=3,cex=0.7)
    lines(sp1[i,j,"inds",],type='b',col=cols1[1])
    lines(2:t,sp2[i,j,"inds",2:t],type='b',col=cols2[1])
    lines(sp10[i,j,"inds",],type='b',lty=2,col=cols1[1])
    lines(2:t,sp20[i,j,"inds",2:t],type='b',col=cols2[1],lty=2)
    points(1:t+.2,sp1[i,j,"seeds",],pch=18,col=cols1[2])
    points(1:t+.2,sp2[i,j,"seeds",],pch=18,col=cols2[2])
    points(1:t+.2,sp10[i,j,"seeds",],pch=5,col=cols1[2],cex=0.7)
    points(1:t+.2,sp20[i,j,"seeds",],pch=5,col=cols2[2],cex=0.7)
    
    #legend("topleft", legend=c("sp1","sp1notlit","seedslit","seedsnolit",
#                               "disp_lit","dispNolit","lit"),col=c(1,"darkred","blue","orange",
 #                                                                  "skyblue","tan",3),lty=1,bty="n")
  }
}
title(main=paste(parms,pnames,sep="=",collapse=" "),outer=T,font.main=1,line=-.3)
dev.off()


# some notes
# dispersal is not dependent on whats there at all..should fix this.
  #like in 1,1 amount dispersed out in timestep 2 is the same. 
  #a constant subtraction but should be a fraction
    #since part of the push mean more seeds -> more disp out 
      #next day (8/6/25) - think I fixed this with difdis function
      #check:
sum(sp1[,,"seeds",2])-sum(sp1[,,"dispered",2])
sum(sp10[,,"seeds",2])-sum(sp10[,,"dispered",2])
sp2[,,"seeds",3]-sp2[,,"dispered",3]
sp20[,,"seeds",3]-sp20[,,"dispered",3]
  #change litter too to use difdis function
      #close but different - nice
#what would bigger rho or bigger l1 look like?
