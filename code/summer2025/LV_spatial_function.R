fig_loc <- "./figures/"
if (!dir.exists(fig_loc)){dir.create(fig_loc)}

source("./dispersal_functions.R")


##function###
runsim <- function(N1, N2, timesteps, r1, r2, alpha11, alpha22, alpha12, alpha21, disp_rate){
	
spparmsvec <- c(r1=r1,r2=r2,alpha11=alpha11,alpha22=alpha22,alpha12=alpha12,alpha21=alpha21)
# Initialize populations
#set.seed(4242)
#runif(grid_size^2, 5, 10)
#N1 <- matrix(rbeta(grid_size^2,0.2,0.2)*10, nrow = grid_size)
#set.seed(4243)
#N2 <- matrix(rbeta(grid_size^2,0.2,0.2)*10, nrow = grid_size)

grid_size <- nrow(N1)

# To store population over time (optional)
pop1_array <- array(NA, dim = c(grid_size, grid_size, timesteps))
pop2_array <- array(NA, dim = c(grid_size, grid_size, timesteps))
growth1 <- pop1_array
growth2 <- pop2_array
pop1_array[,,1] <- N1
pop2_array[,,1] <- N2


# Simulation loop
for (t in 2:timesteps) {
  # Apply Lotka-Volterra competition dynamics
  growth1[,,t] <- r1 * N1 * (1 - alpha11 * N1 - alpha12 * N2)
  growth2[,,t] <- r2 * N2 * (1 - alpha21 * N1 - alpha22 * N2)
  
  N1 <- N1 + growth1[,,t]
  N2 <- N2 + growth2[,,t]
  
  # Ensure populations stay non-negative
  N1[N1 < 0] <- 0
  N2[N2 < 0] <- 0
  
  # Dispersal
  N1 <- disperse2(N1, local_disp=disp_rate, global_disp=0.001)
  N2 <- disperse2(N2, local_disp=disp_rate, global_disp=0.001)
  
  # Store for plotting or analysis
  pop1_array[,,t] <- N1
  pop2_array[,,t] <- N2
  
  #image(N1, main = paste0("Species 1",t), col = terrain.colors(100))
  #image(N2, main = "Species 2", col = terrain.colors(100))
}

#equil reached?##
#equil1 <- mean(pop1_array[,,t]/pop1_array[,,t-1])
equil1 <- max(abs(pop1_array[,,t]/pop1_array[,,t-1]-1))
#equil2 <- mean(pop2_array[,,t]/pop2_array[,,t-1])
equil2 <- max(abs(pop2_array[,,t]/pop2_array[,,t-1]-1))
if(equil1<1e-2){
	cat("pop1 equil reached, ")} else {cat("pop1 non-equil:",equil1,", ")}
if(equil2<1e-2){
	cat("pop2 equil reached")} else {cat("pop2 non-equil:",equil2)}


maxi <- max(c(max(pop1_array),max(pop2_array)))
maxi <- round(maxi) + (5-(round(maxi)%%5))

parms <- paste("S",grid_size,"D",disp_rate,sep="_")
spparms <- paste(names(spparmsvec),spparmsvec,sep="=",collapse='_')
allparms <- paste(parms,spparms)

#fig <- paste0(fig_loc,"popMatD",disp_rate,"S",grid_size,".pdf")
fig <- paste0(fig_loc,"popMatD",allparms,".pdf")
pdf(fig)
par(mfrow=c(grid_size,grid_size),mar=c(1,1,1,1),mgp=c(1,0.5,0),oma=c(3,1,0,0),fg='gray30')
for(i in 1:grid_size){
	for (j in 1:grid_size){
		plot(pop1_array[i,j,], type='l', ylim=c(0,maxi))
		lines(pop2_array[i,j,],col=2, ylim=c(0,maxi))
		if(pop1_array[i,j,1]<pop2_array[i,j,1]){box(col='red3',lty=2)}
		title(main=paste(i,j,sep='-'), line=-0.5, font.main=1, cex.main=0.5)
	}
}
title(sub=allparms,outer=T,line=1)
par(mfrow=c(1,1))
plot(apply(pop1_array,3,sum),type='l',ylim=c(0,maxi*grid_size^2))
lines(apply(pop2_array,3,sum),col=2)
copr <- sum(N1>1 & N2>1)/grid_size^2
dev.off()
return(list(sp1=pop1_array,sp2=pop2_array,copr=copr))
}

#for(i in 1:grid_size){
#	for (j in 1:grid_size){
#plot(growth1[i,j,],type='l')
#lines(growth2[i,j,],col=2)
#}}
#sum(pop$sp1[,,timesteps])
#sum(pop$sp2[,,timesteps])
#pop$sp1[,,1]<pop$sp2[,,1]

# Example plot at final timestep
#par(mfrow = c(1, 2))
#image(N1, main = "Species 1", col = terrain.colors(100))
#image(N2, main = "Species 2", col = terrain.colors(100))


#test dispersal function
#test_mat <- matrix(0.5, nrow = grid_size, ncol = grid_size)
#test_out1 <- disperse(test_mat, disp_rate)
#test_out2 <- disperse(test_mat, disp_rate)
#identical(test_out1, test_out2)  # should be TRUE
#N1 <- matrix(rbeta(grid_size^2,0.2,0.2), nrow = grid_size)
#N2<- N1
#test_out1 <- disperse(N1, disp_rate)
#test_out2 <- disperse(N2, disp_rate)
#identical(test_out1, test_out2) 
#all.equal(test_out1, test_out2) 
