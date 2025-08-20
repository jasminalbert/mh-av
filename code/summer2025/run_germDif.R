#setwd("~/Documents/PhD_school/mh-av/code/summer2025")
timesteps=25
gridsize <- 2
parms <- c("lambda1","lambda2","alpha","beta","rho","theta","l1","l2","D","LLfun","d2")
pnames <- c(lambda10,lambda2,alpha,beta,rho,theta,l1,l2,D,LLfun,decay)
l1 <- 0.9;l2<-0.1
decay <- 0
rho <- 1;theta <- 1
D0 <- 0; D1 <- 0.5
lammax <- 5
lambda10 <- 2
lambda2 <- 2
alpha <- 1.2
beta <- 0.5
beta12 <- alpha; beta21=beta
#LLfun <- "exp" 
LLfun <- "logis"

#individuals at time 1
#N1 <- matrix(c(0.8,0.05,0.1,0.7),ncol=gridsize) 
N10 <- matrix(c(0.8,0,0,0.7),ncol=gridsize) 
N20 <- N10
N20[N10==0] <- rev(N10[N10>0]);N20[N10>0] <-0
#individuals introduced at time tD
tD <- 5
#tD <- 1
N1 <- N20/4; N2 <- N10/4

plot_dispersed <- F
plot_equil <- T
plot_zngi <-F
#fig <- paste0("../../figures/burnin5_lit",decay,".pdf" )
fig <- paste0("../../figures/test3.pdf" )
pdf(fig,width=12)
source("./makesim_germDif.R")
dev.off()

