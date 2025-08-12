timesteps=10
gridsize <- 2
parms <- c("lambda1","lambda2","alpha","beta","rho","theta","l1","l2","D")
l1 <- 0.9;l2<-0.1
rho <- 0.9;theta <- 1
D <- 0.5
lambda1 <- 2
lambda2 <- 2
alpha <- 0.5
beta <- 0.5
#make a loop
#individuals at time 1
N1 <- matrix(c(0.8,0.05,0.1,0.7),ncol=gridsize) 
N2 <- 0.85-N1

source("./makesim_germDif.R")