
#litter differences
litDif <- function(litter){
  gridsize <- dim(litter)[1]
  zeropad <- matrix(0,nrow=gridsize+2,ncol=gridsize+2)
  zeropad[(1:gridsize)+1,(1:gridsize)+1] <- litter
  litdif <- ifelse(is.na(litter),NA,0)
  for (i in 1:gridsize){
    for (j in 1:gridsize){
      nbs <- expand.grid(i=(i-1):(i+1),j=(j-1):(j+1))
      nbs <- nbs[!(nbs$i==i&nbs$j==j),]
      for (k in 1:nrow(nbs)){
        #nb <- litter[nbs$i[k],nbs$j[k]]#;print(nb) 
        nb <- zeropad[nbs$i[k]+1,nbs$j[k]+1]#;print(nb) 
        if (nbs$i[k]%in%(1:nrow(litter)) & nbs$j[k]%in%(1:ncol(litter))){
          dif <- nb-litter[i,j]
        } else {dif <- 0}
        #dif <- ifelse(length(nb)>0,nb-litter[i,j],0)
        litdif[i,j] <- litdif[i,j]+dif
        #cat("nb",nb,"dif",dif,"litdif",litdif[i,j],"\n")
      }
    }
  }
  return(litdif)
}
litDif(litter)
litdif
litter
(0.05-0.4)+(-0.4+0.35)+(-0.4+0.025) #[1,1]: -.775
(0.4-0.025)+(0.05-0.025)+(0.35-0.025) #[2,1]: .725
(0.4-0.05)+(0.025-0.05)+(0.35-0.05) #[1,2]: .625
(0.4-0.35)+(0.025-0.35)+(0.05-0.35) #[2,2]: -.575
#kernel <- matrix(1,nrow=3,ncol=3)
#litdif[i,j] <- litdif[i+1,j]
#i=2;j=1
#litdif[(i-1):(i+1),(j-1):(j+1)]  
#c(litdif[nbs$i,nbs$j]  )
#is.integer(litdif[i-1,j])





gridsize <- 2
N1 <- matrix(c(0.8,0.05,0.1,0.7),ncol=gridsize) #individuals at time 1

sp1array <- array(NA,dim=c(gridsize,gridsize,timesteps))

sp1array[,,1] <- N1; #individuals at time 1
N1seeds <- lambda1*N1*(1-N1) #individuals make seeds
#seeds disperse
dmat <- matrix(rnorm(gridsize^2),ncol=gridsize)
N1seeds_dispersed <- N1seeds*0.7+sum(N1seeds*0.3)*dmat/sum(dmat) #come back to fill this out
sum(N1seeds_dispersed)-sum(N1seeds)<1e-6
lit1 <- N1*0.5 #adults die and turn into litter. litter made in time 1 impacts seeds in time 2
litter <- lit1+0
litdif <- litDif(litter)
#new dispersed seeds represent that germinants->adults of time 2
N1_2 <- N1seeds_dispersed
#again they compete and reproduce, but now there is litter from tme 1
#litdif*lambda1 #litdif<0 means other cells have less litter than you
#for mh, its a push, positive feedback = reverse diffusion
rho<-0.5
lambda1lit <- lambda1*(1-rho*litdif)
N1_2seeds <- lambda1*N1_2*(1-N1_2)
N1_2seedslit <- lambda1lit*N1_2*(1-N1_2)

easydisperse <- function(N1seeds,seed,D=0.3){
  set.seed(seed)
  dmat <- matrix(rnorm(gridsize^2),ncol=gridsize)
  N1seeds_dispersed <- N1seeds*(1-D)+sum(N1seeds*D)*dmat/sum(dmat) #come back to fill this out
  print(sum(N1seeds_dispersed)-sum(N1seeds)<1e-6)
  return(N1seeds_dispersed)
}

#make a loop
timesteps=10
gridsize <- 2
N1 <- matrix(c(0.8,0.05,0.1,0.7),ncol=gridsize) #individuals at time 1
sp1array <- array(NA,dim=c(gridsize,gridsize,timesteps))
seedsnolit <- sp1array
dispersednolit <- sp1array
dispersedlit <- sp1array
seedslit <- sp1array
lit <- sp1array
lit[,,1] <- ifelse(is.na(N1),NA,0)
sp1array[,,1] <- N1; #individuals at time 1
sp1arraynolit <- sp1array
t=1
for (t in 1:timesteps){
  seed=sample(1:1e6,1)
  litdif <- litDif(lit[,,t])
  lambda1lit <- lambda1*(1-rho*litdif)
  seedslit[,,t] <- lambda1lit*sp1array[,,t]*(1-sp1array[,,t])
  seedsnolit[,,t] <- lambda1*sp1arraynolit[,,t]*(1-sp1arraynolit[,,t])
  dispersedlit[,,t] <- easydisperse(seedslit[,,t],seed=seed)  
  dispersednolit[,,t] <- easydisperse(seedsnolit[,,t],seed=seed)  
  lit[,,t+1] <- sp1array[,,t]*0.5
  sp1array[,,t+1] <- dispersedlit[,,t]
  sp1arraynolit[,,t+1] <- dispersednolit[,,t]
}
par(mfrow=c(2,2),xpd=F)
timesteps=10;t=10
for (i in 1:2){
  for (j in 1:2){
    plot(1:timesteps,seq(0,1,length.out=t),type="n",ylab='',xlab='',xlim=c(1,7))
    lines(1:timesteps,sp1array[i,j,],type='b')
    lines(1:timesteps,sp1arraynolit[i,j,],type='b',col='darkred')
    lines(1:timesteps+.3,seedsnolit[i,j,],type='b',col='orange')
    lines(1:timesteps+.3,seedslit[i,j,],type='b',col='blue')
    #lines(1:timesteps+.5,dispersedlit[i,j,],type='b',col='skyblue',lty=2)
    #lines(1:timesteps+.5,dispersednolit[i,j,],type='b',col='tan',lty=2)
    lines(1:timesteps-.2,lit[i,j,],col=3,type='b')
    #points(2,N1_2[i,j])
    #points(2.3,N1_2seeds[i,j])
    #points(2.3,N1_2seedslit[i,j],pch=5)
    legend("topleft", legend=c("sp1","sp1notlit","seedslit","seedsnolit",
          "disp_lit","dispNolit","lit"),col=c(1,"darkred","blue","orange",
              "skyblue","tan",3),lty=1,bty="n")
  }
}
# try + lifdif^2?
i=1;j=1



#species 2
#make a loop
timesteps=5
gridsize <- 2
N2 <- 0.85-N1
sp2array <- array(NA,dim=c(gridsize,gridsize,timesteps))
seedsnolit <- sp2array
dispersednolit <- sp2array
dispersedlit <- sp2array
seedslit <- sp2array
#lit <- sp1array
#lit[,,1] <- ifelse(is.na(N1),NA,0)
sp2array[,,2] <- N2; #individuals at time 1
sp2arraynolit <- sp2array
D=0.3;t=2
for (t in 2:timesteps){
  seed=sample(1:1e6,1)
  litdif <- litDif(lit[,,t])
  litdifnorm <- litdif/sum(abs(litdif));
    #lambda1lit <- lambda1*(1-rho*litdif)
  seedslit[,,t] <- lambda2*sp2array[,,t]*(1-sp2array[,,t])
  seedsnolit[,,t] <- lambda2*sp2arraynolit[,,t]*(1-sp2arraynolit[,,t])
  dispersedlit[,,t] <- seedslit[,,t] - D*litdifnorm
  dispersednolit[,,t] <- easydisperse(seedsnolit[,,t],seed=seed)  
  #lit[,,t+1] <- sp1array[,,t]*0.5
  sp2array[,,t+1] <- dispersedlit[,,t]
  sp2arraynolit[,,t+1] <- dispersednolit[,,t]
}
t=1+t
#lif dif = sum(neighbors-you)
#positive: nbs have more than you, you disperse out(reverse diffusion)
#negative: you have more than neighbors, they disp in (rev diff)
par(mfrow=c(1,4))
for(i in 1:2){for(j in 1:2){plot(5:10,lit[i,j,5:10],main=paste(i,j),type='b',ylim=c(-0.65,0.65))
  points(10,litdif[i,j],col=3,pch=19);abline(h=c(0,lit[i,j,10]),lty=c(1,3))
  points(10,litdifnorm[i,j],col='magenta')}}
((seedslit*litdifnorm)/sum(seedslit*litdifnorm))*seedslit
#litdifnorm : movement of your cell
D=0.3
par(mfrow=c(2,2),xpd=F)
timesteps=5;t=5
for (i in 1:2){
  for (j in 1:2){
    plot(1:timesteps,seq(0,1,length.out=t),type="n",ylab='',xlab='',xlim=c(1,t))
    lines(1:timesteps,sp2array[i,j,],type='b')
    lines(1:timesteps,sp2arraynolit[i,j,],type='b',col='darkred')
    lines(1:timesteps+.3,seedsnolit[i,j,],type='b',col='orange')
    lines(1:timesteps+.3,seedslit[i,j,],type='b',col='blue')
    lines(1:timesteps+.5,dispersedlit[i,j,],type='b',col='skyblue',lty=2)
    lines(1:timesteps+.5,dispersednolit[i,j,],type='b',col='tan',lty=2)
    lines(1:timesteps-.2,lit[i,j,1:timesteps],col=3,type='b')
    legend("topleft", legend=c("sp1","sp1notlit","seedslit","seedsnolit",
                               "disp_lit","dispNolit","lit"),col=c(1,"darkred","blue","orange",
                                                                   "skyblue","tan",3),lty=1,bty="n")
  }
}
plot(1:timesteps,seq(0,1,length.out=t),type="n",ylab='',xlab='')
lines(1:timesteps,sp1array[i,j,])
lines(1:timesteps+.3,seeds[i,j,])
lines(1:timesteps+.3,seedslit[i,j,],lty=2, type='b')
lines(1:timesteps+.5,dispersed[i,j,],col="lightblue", lty=3,type='b')
lines(1:timesteps+.8,lit[i,j,],col=3)
par(mfrow=c(2,2))
t=3
for (i in 1:2){
  for (j in 1:2){
    plot(1:t,seq(0,1,length.out=t),type="n",ylab='',xlab='')
    points(1,N1[i,j])
    points(1.3,N1seeds[i,j])
    points(1.5,N1seeds_dispersed[i,j])
    points(1.8,litter[i,j],col=3)
    points(2,N1_2[i,j])
    points(2.3,N1_2seeds[i,j])
    points(2.3,N1_2seedslit[i,j],pch=5)
  }
}

litdifnorm
sum(litdifnorm);sum(abs(litdifnorm))
disp <-;sum(disp);sum(seedslit)
disp; seedslit

litdifnorm <- litdif/sum(abs(litdif))
disp <- seedslit*(1-litdifnorm) + seedslit*litdifnorm
sum(disp)
sum(seedslit) #the same right now

