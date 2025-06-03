source("./litter_functions.r")

decay_half_life <- c(3,0.5)
decay_rate <- log(2) / decay_half_life
timesteps <- 30
# Initialize
litter_cohorts <- initialize_litter_cohorts(grid_size)

lit_array <- array(NA, dim = c(grid_size, grid_size, timesteps))
for (t in 1:timesteps) {
  # Compute new litter (e.g., based on N1 (mh), N2(av) losses)
  #new_litter <- 0.1 * N1 + 0.3 * N2  # example
  N1 <- matrix(rbeta(grid_size^2,0.2,0.2)*10, nrow = grid_size)
  N2 <- matrix(rbeta(grid_size^2,0.2,0.2)*10, nrow = grid_size)

  # Add new litter
  litter_cohorts <- add_litter(litter_cohorts, N1,N2, t)

  # Update and get current total litter
  decay_result <- decay_litter(litter_cohorts, t, decay_rate[1], decay_rate[2],grid_size)
  total_litter <- decay_result$litter
  litter_cohorts <- decay_result$updated_cohorts

  # Use total_litter in your model, if needed
  lit_array[,,t] <- total_litter
}
#lapply(litter_cohorts, colSums)

par(mfrow=c(grid_size,grid_size),mar=c(1,1,1,1),mgp=c(1,0.5,0),oma=c(3,1,0,0),fg='gray30')
for(i in 1:grid_size){
	for (j in 1:grid_size){
		plot(lit_array[i,j,], type='l',ylim=c(0,7))
		lit <- litter_cohorts[[ij_to_index(i,j,grid_size)]]
		points(lit$birth_time,lit$amount1,pch=19,cex=.4)
		points(lit$birth_time,lit$amount2,col=2,pch=19,cex=.4)
		#cat(i,j,1,cor(lit_array[i,j,],litter_cohorts[[ijx]]$amount1))
		#cat(2,cor(lit_array[i,j,],litter_cohorts[[ijx]]$amount2),"\n")
	}
}


