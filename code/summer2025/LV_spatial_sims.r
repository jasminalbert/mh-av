source("./LV_spatial_function.R")
# Species parameters
r1 <- 0.1
r2 <- 0.1
alpha11 <- 0.06
alpha22 <- 0.06
alpha12 <- 0.08
alpha21 <- 0.08
# Dispersal rate (0 = no dispersal, 1 = full mixing)
disp_rate <- 0.1

# SIM Parameters
grid_size <- 5 #S
timesteps <- 300

S_vec <- c(5,10,15,20,25) #grid sizes
copr <- rep(NA,length(S_vec))
maxS <- max(S_vec)
set.seed(33)
N1pool <- matrix(rbeta(maxS^2,0.2,0.2)*10, nrow = maxS)
N2pool <- matrix(rbeta(maxS^2,0.2,0.2)*10, nrow = maxS)

for (s in seq_along(S_vec)){
	grid_size <- S_vec[s]
	N1 <- N1pool[1:grid_size,1:grid_size]
	N2 <- N2pool[1:grid_size,1:grid_size]
	cat("S:",grid_size," ")
	#set.seed(1234)
	pop<-runsim(N1,N2, timesteps, r1, r2, alpha11, alpha22, alpha12, alpha21, disp_rate)
	copr[s] <-pop$copr
	cat("\n")
	print(sum(pop$sp1[,,timesteps])/(grid_size^2))
	print(sum(pop$sp2[,,timesteps])/(grid_size^2))
};print(copr);
#confused why winner is the same


#need to add litter feed backs
