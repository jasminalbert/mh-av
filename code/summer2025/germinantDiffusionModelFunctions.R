
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








