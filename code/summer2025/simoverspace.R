setwd("/Users/jasminalbert/Documents/PhD_school/mh-av/code/summer2025")
source("./dispersal_functions.R")
source("./litter_functions.r")
source("./LV_space_wLitter.R")
coexprob <- function(sp1,sp2){
	Ntot <- sp1+sp2
	prob1 <- sp1/Ntot
	prob2 <- sp2/Ntot
	return(sum(prob1<0.8 & prob2<0.8)/length(Ntot))
}
plotcoprob <- function(poplist, new=T){
	gridsizes <- as.numeric(names(poplist))
	cprob <- {}
	if (new){
		par(mfrow=c(1,1),mar=c(1,1,1,1), mgp=c(1.2,0.5,0),oma=c(3,1,0,0),fg='gray30', oma=c(2,8,1,8),xpd=NA)
		plot(gridsizes,ylim=c(0,1),type="n",ylab='coexistence probability', xlab="grid_size")
		
		parms <- poplist[[1]]$parms
		pop.ptxt <- paste(names(parms[[1]]),parms[[1]],sep="=",collapse=" ")
		lit.ptxt <- paste(names(parms[[2]]),parms[[2]],sep="=",collapse=" ")
		disp.ptxt <- paste(names(parms[[3]]),parms[[3]],sep="=",collapse=" ")
		text(10,1,labels=pop.ptxt,adj=1)
		text(10,0.975,labels=lit.ptxt,adj=1)
		text(10,0.95,labels=disp.ptxt,adj=1)
	}
	for(grid in gridsizes){
		grid_size <- grid;pop <- poplist[[grid]]
		#coexistence probability
		sp1 <- pop$sp1[,,timesteps];sp2 <- pop$sp2[,,timesteps]
		cprob[grid] <- coexprob(sp1,sp2)
		points(grid_size,cprob[grid],cex=2)
	}
	lines(gridsizes,cprob,col=rgb(0,0,0,0.4))
	return(cprob)
}

plotpatches <- function(poplist,ylim=default){
	gridsizes <- as.numeric(names(poplist))
	for(grid in gridsizes){
		grid_size <- grid
		pop <- poplist[[as.character(grid)]]
		par(mfrow=c(grid_size,grid_size),mar=c(1,1,1,1), mgp=c(1,0.5,0),oma=c(3,1,0,0),fg='gray30')
		maxi <- max(c(max(pop$sp1),max(pop$sp2)),na.rm=T)
		maxi <- round(maxi) + (5-(round(maxi)%%5))
#maxi = 100
		ymax <- ylim
		if(ylim=="default"){ymax <- maxi}
		for(i in 1:grid_size){
			for (j in 1:grid_size){
				plot(pop$sp1[i,j,], type='l', ylim=c(0,ymax),xlim=c(0,timesteps),xlab='',ylab='')
				lines(pop$sp2[i,j,],col=2)
				lines(pop$litter[i,j,],col="lightgrey")
				if(pop$sp1[i,j,1]<pop$sp2[i,j,1]){box(col='red3',lty=2)}
				title(main=paste(i,j,sep='-'), line=-0.5, font.main=1, cex.main=0.5)
				pop1d <- pop$sp1[i,j,c(1,timesteps)]
				pop2d <- pop$sp2[i,j,c(1,timesteps)]
				text(rep(c(.15,1),2)*timesteps,y=c(pop1d,pop2d)+maxi*.07, labels=round(c(pop1d,pop2d),2),adj=1,col=rep(1:2,each=2), cex=0.7)
			}
		}	
	}
}


nsims <- 100
nsims <- 20
simslist <- list()
gridsizes <- 1:10
gridsizes=5
gridnames <- as.character(gridsizes)
cprobmat <- array(dim=c(nsims,length(gridsizes)))
set.seed(3336)
seepatches <- T
for (sim in 1:nsims){
	poplist <- vector(mode="list",length=length(gridsizes))
	names(poplist) <- gridnames 
	#gridsizes <- as.numeric(names(poplist))
	
	N2all <- rbeta(max(gridsizes)^2,0.2,0.2)*10
	N1all <- rbeta(max(gridsizes)^2,0.2,0.2)*10
	
	for (grid in gridnames){
		grid_size <- as.numeric(grid)
	
		#set.seed(3336)
		N2 <- matrix(N2all[1:grid_size^2], nrow=grid_size)
		N1 <- matrix(N1all[1:grid_size^2], nrow=grid_size)
		
		poplist[[grid]] <- runsim(timesteps, N2,N1,popparms,dispparms,litparms)
	}
	if (seepatches){plotpatches(poplist[grid])}
	poplist[[1]]$parms
	simslist[[sim]] <- poplist
	cprobmat[sim,]<-plotcoprob(poplist, new=ifelse(sim==1,T,F))
	cat("\nsim",sim,"done...")
	#if (sim){
	#	plotcoprob(poplist)
	#} else {plot}
}
poplist[[1]]$parms
saveRDS(simslist,"simslist_neutral.RDS")
saveRDS(cprobmat,"coexprob_neutral.RDS")
cprobmat<-readRDS("coexprob_neutral.RDS")
simslist<-readRDS("simslist_neutral.RDS")
dim(cprobmat)
colMeans(cprobmat)
cprobmat <- cprobmat[1:20,]
simslist[[20]][[10]]$sp2[10,10,50]+simslist[[20]][[10]]$sp1[10,10,50]
simslist[[20]][[10]]$sp2[10,10,50]/0.5297117
pop <- simslist[[1]][[1]]
pdf("cprobfigneutral_20sims.pdf")
txtr <- 15
txtu <- 1.06
for (sim in 1:nsims){
	if (sim==1){
		par(mfrow=c(1,1),mar=c(1,1,1,1), mgp=c(1.2,0.5,0),oma=c(3,1,0,0),fg='gray30', oma=c(2,8,1,8),xpd=NA)
		plot(gridsizes,ylim=c(0,1),type="n",ylab='coexistence probability', xlab="grid_size",bty="l")
		#pop <- poplist[[1]]
		parms <- pop$parms
		pop.ptxt <- paste(names(parms[[1]]),parms[[1]],sep="=",collapse=" ")
		lit.ptxt <- paste(names(parms[[2]]),parms[[2]],sep="=",collapse=" ")
		disp.ptxt <- paste(names(parms[[3]]),parms[[3]],sep="=",collapse=" ")
		text(txtr,txtu,labels=pop.ptxt,adj=1)
		text(txtr,txtu-0.03,labels=lit.ptxt,adj=1)
		text(txtr,txtu-0.06,labels=disp.ptxt,adj=1)
	}
	#for(grid in gridsizes){
	#	grid_size <- grid;pop <- poplist[[grid]]
		#coexistence probability
	#	sp1 <- pop$sp1[,,timesteps];sp2 <- pop$sp2[,,timesteps]
	#	cprob[grid] <- coexprob(sp1,sp2)
	#	points(grid_size,cprob[grid],cex=2)
	#}
	lines(gridsizes,cprobmat[sim,],col=rgb(0,0,0,0.2))
}
lines(gridsizes,colMeans(cprobmat),col=3,lwd=3)	
dev.off()

set.seed(3336)
seepatches <- T
poplist <- vector(mode="list",length=length(gridsizes))
names(poplist) <- gridnames 
	#gridsizes <- as.numeric(names(poplist))
	
N2all <- rbeta(max(gridsizes)^2,0.2,0.2)*10
N1all <- rbeta(max(gridsizes)^2,0.2,0.2)*10
	

grid_size <- as.numeric(grid)
	
		#set.seed(3336)
N2 <- matrix(N2all[1:grid_size^2], nrow=grid_size)
N1 <- matrix(N1all[1:grid_size^2], nrow=grid_size)
timesteps=20
poplist[[grid]] <- runsim(timesteps, N2,N1,popparms,dispparms,litparms)
if (seepatches){plotpatches(poplist[grid], ylim=100)}
poplist[[1]]$parms

pdf("finallyPE.pdf")
plotpatches(poplist[grid])
dev.off()

#n <- 1
#grid <- gridsizes[n];n=n+1