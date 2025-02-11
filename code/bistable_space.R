# Load required packages
library(ggplot2)

# Parameters
grid_size <- 50  # Size of the spatial grid (50x50)
timesteps <- 100 # Number of simulation steps
switch_prob <- 0.1  # Base probability of switching states
heterogeneity_strength <- 0.3  # Influence of spatial heterogeneity

# Initialize grid with bistable states (A or B)
set.seed(42)
grid <- matrix(sample(c("A", "B"), grid_size^2, replace = TRUE), nrow = grid_size)

# Generate spatial heterogeneity (random field)
heterogeneity <- matrix(runif(grid_size^2, -1, 1), nrow = grid_size)

# Function to update the grid based on local interactions
update_grid <- function(grid, heterogeneity, switch_prob) {
  new_grid <- grid  # Copy current grid
  
  for (i in 1:nrow(grid)) {
    for (j in 1:ncol(grid)) {
      # Count neighbors in state "A"
      neighbors <- c(
        grid[max(i-1, 1), j], grid[min(i+1, nrow(grid)), j],
        grid[i, max(j-1, 1)], grid[i, min(j+1, ncol(grid))]
      )
      num_A <- sum(neighbors == "A")
      num_B <- sum(neighbors == "B")
      
      # Environmental influence
      env_factor <- heterogeneity[i, j] * heterogeneity_strength
      #make heterogeneity litter
      #and species specific HS (influence of litter on rec)
      
      #base switch prob should be low, they dont switch unless theres litter
      #or is it outcompete prob?
      #av can outcompete mh in high litter so its still litter
      
      
      # Probability of switching states
      # B is favored by env 
      # env f influences switch from A to B
      # and help keep B in place
      # so B is avena
      # having neighbors with same state helps maintain state?
      # probably true bc dispersal.
      prob_A_to_B <- switch_prob + env_factor - (num_A / 4) * 0.2
      prob_B_to_A <- switch_prob - env_factor - (num_B / 4) * 0.2
      # Probability of staying the same is runif(1)?
      
      # Apply state transitions
      if (grid[i, j] == "A" && runif(1) < prob_A_to_B) {
        new_grid[i, j] <- "B"
      } else if (grid[i, j] == "B" && runif(1) < prob_B_to_A) {
        new_grid[i, j] <- "A"
      }
    }
  }
  return(new_grid)
}

# Run the simulation
for (t in 1:timesteps) {
  grid <- update_grid(grid, heterogeneity, switch_prob)
}

# Convert grid to data frame for visualization
grid_df <- expand.grid(x = 1:grid_size, y = 1:grid_size)
grid_df$state <- as.vector(grid)

# Plot final state of the landscape
ggplot(grid_df, aes(x = x, y = y, fill = state)) +
  geom_tile() +
  scale_fill_manual(values = c("A" = "blue", "B" = "yellow")) +
  theme_minimal() +
  labs(title = "Bistable Landscape Simulation", x = "X Coordinate", y = "Y Coordinate")
  
gridlist <- list()  
grid_df <- expand.grid(x = 1:grid_size, y = 1:grid_size)
grid2 <- array(NA, dim=dim(grid))
# Plot every timestep
switch_prob = 0.5
for (t in 1:timesteps) {
  grid <- update_grid(grid, heterogeneity, switch_prob)
  gridlist[[t]] <- grid
grid2[grid=="A"] <- which(LETTERS=='A')
grid2[grid=="B"] <- which(LETTERS=='B')
image(grid2, main=paste("time=",t))
}  
grid = gridlist[[t]]





