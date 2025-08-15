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

difdis <- function(rd,seeds){
	dispin <- rd*seeds
	dispout=(-rd*seeds)/3
	dis=array(NA,dim=dim(seeds))
	for (i in 1:(ncol(seeds)^2)){
		dis[i] = seeds[i]-sum(dispout[-i]) #whats left after you gave
		#dis_1[i] = seeds[i]-rd[i]*seeds[i]
	}
	res <- dis - dispin
	return(res)
}

#randomdispersalcheck.R for sanity check 


###LOGISTIC FUNCTION#

#L:carry capacity (max lambda)
#k: logistic growth rate, steepness of curve (rho)
#x0: x value at functions midpoint (lambdai^0)
#standard logistic function has L=k=1,x0=0
logis <-function(x,L=1,k=1,x0=0){
  denom <- 1+exp(-k * (x-x0))
  y <- L/denom
  return(y)
}




