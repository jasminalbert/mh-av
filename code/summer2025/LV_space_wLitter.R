source("./dispersal_functions.R")
source("./litter_functions.r")

#args
#timesteps

#N1			#init matrix
#N2			#init matrix
#popparms	#c(r1,r2,alpha11,alpha22,alpha12,alpha21,g1,g2) 
#dispparms	#c(d1,d2,gd1,gd2)
#litparms	#c(l1,l2,hali1,hali2,theta1,theta2,rho1,rho2)
runsim <- function(timesteps, N1, N2, popparms, dispparms, litparms){
	#parameters
	list2env(as.list(popparms), envir = environment())
	r1 <- popparms["r1"]
	r2 <- popparms["r2"]
	alpha11 <- popparms["alpha11"]
	alpha22 <- popparms["alpha22"]
	alpha12 <- popparms["alpha12"]
	alpha21 <- popparms["alpha21"]
	g1 <- popparms["g1"]
	g2 <- popparms["g2"]
	
	list2env(as.list(litparms), envir=environment())
	list2env(as.list(dispparms), envir=environment())
	parms <- list(popparms, litparms, dispparms)
	
	#initialize litter structure
	litter_cohorts <- initialize_litter_cohorts(grid_size)
	decay_rate <- log(2) / c(hali1,hali2)
	
	#store pop
	grid_size <- nrow(N1)
	pop1_array <- array(NA, dim = c(grid_size, grid_size, timesteps))
	pop2_array <- array(NA, dim = c(grid_size, grid_size, timesteps))
	lit_array <- pop1_array
	growth1 <- pop1_array
	growth2 <- pop2_array
	pop1_array[,,1] <- N1
	pop2_array[,,1] <- N2
	
	# First litter
	litter_cohorts <- add_litter(litter_cohorts, N1,N2, time_step=1)
	#print(litter_cohorts)
	decay_result <- decay_litter(litter_cohorts, current_time=1, decay_rate[1], 	decay_rate[2],grid_size) #need these lines here?
    total_litter <- decay_result$litter
    litter_cohorts <- decay_result$updated_cohorts
    lit_array[,,1] <- total_litter
	
	#LV pop sim
	for (time in 2:timesteps) {
		cat(time)
		# Dispersal: seeds from previous time step disperse
  		N1 <- disperse2(N1, local_disp=d1, global_disp=gd1)
  		N2 <- disperse2(N2, local_disp=d2, global_disp=gd2)
  		
  		# Litter responses 
  		r10 <- min((total_litter*rho1+1)*r1,r2*1.5)
  		r20 <- min((total_litter*rho2+1)*r2,r1*1.5)
  		g10 <- min((total_litter*theta1+1)*g1,1)
  		g20 <- min((total_litter*theta2+1)*g2,1)
  		
  		# Apply Lotka-Volterra competition dynamics
  		growth1[,,time] <- g10 * r10 * N1 * (1 - alpha11 * g10 * N1 - alpha12 * g20 * N2)
  		growth2[,,time] <- g20 * r20 * N2 * (1 - alpha21 * g10 * N1 - alpha22 * g20 * N2)
  
  		N1 <- growth1[,,time] #+N1
  		N2 <- growth2[,,time] #+N2
		#print(growth1[,,time])
		#print(N1);print(N2)  
  		# Ensure populations stay non-negative
  		N1[N1 < 0] <- 0
  		N2[N2 < 0] <- 0
  
  		
  		# Track litter
  		litter_cohorts <- add_litter(litter_cohorts, N1,N2, time_step=time)
		decay_result <- decay_litter(litter_cohorts, current_time=time, decay_rate[1], 	decay_rate[2],grid_size) 
    	total_litter <- decay_result$litter
    	litter_cohorts <- decay_result$updated_cohorts
    	
    	
  		# Store for plotting or analysis
  		lit_array[,,time] <- total_litter
  		pop1_array[,,time] <- N1
  		pop2_array[,,time] <- N2
  
  		#image(N1, main = paste0("Species 1",t), col = terrain.colors(100))
  		#image(N2, main = "Species 2", col = terrain.colors(100))
	}
return(list(sp1=pop1_array,sp2=pop2_array,litter=lit_array,litco=litter_cohorts,parms=parms))

}
#pop <- list(sp1=pop1_array,sp2=pop2_array,litter=lit_array,litco=litter_cohorts)

r1 <- 0.25
r2 <- 0.2# 0.05
alpha11<- 0.08#0.02
alpha22<- 0.08#2
alpha12<- 0.25
alpha21<- 0.2#6
g1<- 0.5
g2<- 0.2#1
popparms <- c(r1=r1,r2=r2,alpha11=alpha11,alpha22=alpha22, alpha12=alpha12,alpha21=alpha21,g1=g1,g2=g2)
popparms <- c(r1=r1,r2=r2,alpha11=alpha11,alpha22=alpha22, alpha12=0.06,alpha21=0.06,g1=g1,g2=g2)

l1 <- 0.1
l2 <- 0.1#2 #AV's N to litter ratio
hali1 <- 4#4 #half life
hali2 <- 0.5
theta1 <- 0
theta2 <- 0.6#0.6 #AV's germination response to litter
rho1 <- 0.1
rho2 <- 0 #AV's fecundity response to litter
litparms <-c(l1=l1,l2=l2,hali1=hali1,hali2=hali2,theta1=theta1,theta2=theta2,rho1=rho1,rho2=rho2)
litparms <-c(l1=l1,l2=l2,hali1=hali1,hali2=hali2,theta1=0,theta2=0,rho1=0,rho2=0)

d1 <- 0.1
d2 <- 0.1
gd1 <- 0.001
gd2 <- 0.001
dispparms <- c(d1=d1,d2=d2,gd1=gd1,gd2=gd2)
dispparms <- c(d1=0,d2=0,gd1=0,gd2=0)
dispparms <- c(d1=0.0075,d2=0.0075,gd1=0.001,gd2=0.001)
grid_size <- 5
timesteps <- 50

#set.seed(333)
#N2 <- matrix(rbeta(grid_size^2,0.2,0.2)*10, nrow = grid_size)
#N1 <- matrix(rbeta(grid_size^2,0.2,0.2)*10, nrow = grid_size)

#pop <- runsim(timesteps, N2,N1,popparms,dispparms,litparms)

#par(mfrow=c(grid_size,grid_size),mar=c(1,1,1,1),mgp=c(1,0.5,0),oma=c(3,1,0,0),fg='gray30')
# maxi <- max(c(max(pop$sp1),max(pop$sp2)),na.rm=T)
# maxi <- round(maxi) + (5-(round(maxi)%%5))
# #maxi = 100
# for(i in 1:grid_size){
	# for (j in 1:grid_size){
		# plot(pop$sp1[i,j,], type='l', ylim=c(0,30),xlim=c(0,timesteps))
		# lines(pop$sp2[i,j,],col=2)
		# lines(pop$litter[i,j,],col="lightgrey")
		# if(pop$sp1[i,j,1]<pop$sp2[i,j,1]){box(col='red3',lty=2)}
		# title(main=paste(i,j,sep='-'), line=-0.5, font.main=1, cex.main=0.5)
	# }
# }
# #pop$litter[,,timesteps]