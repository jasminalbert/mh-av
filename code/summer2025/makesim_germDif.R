#working on (1) burn in monoculture

#storage array
dimnames <- list(i=1:gridsize, j=1:gridsize,stage=c("inds","seeds","dispered"), 
                 time=1:(timesteps+1))
sp1 <- array(NA, dim=c(gridsize,gridsize,3,timesteps+1), dimnames)
sp2 <- sp1
#litter
lit <- sp1[,,1,]
#saving lambdalit
lambda <- lit
#initializing
sp1[,,"inds",1] <- N10
sp2[,,"inds",1] <- N20
#sp2[,,"inds",2] <- N2
lit[,,1] <- ifelse(is.na(N1),NA,0)

#no litter models
sp10 <- sp1; sp20 <- sp2
t=1
for (t in 1:timesteps){
  if (t<tD){D<-D0} else{D<-D1}
  if (t==tD){
    sp1[,,"inds",t] <- N1+sp1[,,"inds",t]
    sp10[,,"inds",t] <- N1+sp10[,,"inds",t]
    sp2[,,"inds",t] <- N2+sp2[,,"inds",t]
    sp20[,,"inds",t] <- N2+sp20[,,"inds",t]
  }
  #no lit at t=1
  litdif <- litDif(lit[,,t])
  if (LLfun =="logis"){lambda1lit <- logis(-litdif,L=(lambda10-1)*2,k=rho)+1}
  if (LLfun == "exp"){lambda1lit <- lambda10*exp(-litdif*rho)}
  sp1[,,"seeds",t] <- lambda1lit*sp1[,,"inds",t]*(1-sp1[,,"inds",t]-alpha*sp2[,,"inds",t])
  sp10[,,"seeds",t] <- lambda10*sp10[,,"inds",t]*(1-sp10[,,"inds",t]-alpha*sp20[,,"inds",t])
  sp2[,,"seeds",t] <- lambda2*sp2[,,"inds",t]*(1-sp2[,,"inds",t]-beta*sp1[,,"inds",t])
  sp20[,,"seeds",t] <- lambda2*sp20[,,"inds",t]*(1-sp20[,,"inds",t]-beta*sp10[,,"inds",t])
  #zero <- mapply(function(x){x[,,"seeds",t]},list(sp1,sp10,sp2,sp20))<0
  #if (any(zero)){cat("\nZERO!",t)}
  #mapply(function(x){x[,,"seeds",t]},list(sp1,sp10,sp2,sp20))
  if(any(na.omit(sp1[,,"seeds",t]) <0)){sp1[,,"seeds",t][na.omit(sp1[,,"seeds",t])<0]<-0}
  if(any(na.omit(sp2[,,"seeds",t]) <0)){sp2[,,"seeds",t][na.omit(sp2[,,"seeds",t])<0]<-0}
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
  #if (t>1){
  sp2[,,"inds",t+1] <- sp2[,,"dispered",t]
  sp20[,,"inds",t+1] <- sp20[,,"dispered",t]
  #}
  #save lambda
  lambda[,,t] <- lambda1lit
}
litmaxes <- apply(lit[,,-1], 3, function(x) which(x==max(x),T ))
litmax <- lit
litmax[,,] <- FALSE
litmax[cbind(t(litmaxes),colnames(litmaxes))] <- TRUE

##mapply(function(x){x<0}, list(sp1[,,"dispered",],sp10[,,"dispered",],sp2[,,"dispered",],sp20[,,"dispered",] ))

source("./plotsim.R")
