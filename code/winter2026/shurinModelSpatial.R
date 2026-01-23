get_neighbors <- function(mat, i, j) {
  ni <- nrow(mat)
  nj <- ncol(mat)
  
  rows <- pmax(1, i - 1):pmin(ni, i + 1)
  cols <- pmax(1, j - 1):pmin(nj, j + 1)
  
  mat[rows, cols]
}
simulate_lattice <- function(
    n = 50,
    timesteps = 200,
    ci = 0.4, cj = 0.4,
    ei = 0.1, ej = 0.1,
    ri = 0.5, rj = 0.3
) {
  
  # Initialize lattice (mostly empty)
  grid <- matrix(1, n, n)
  
  # Seed species
  grid[sample(length(grid), n^2 * 0.05)] <- 2
  grid[sample(length(grid), n^2 * 0.05)] <- 3
  
  history <- vector("list", timesteps)
  
  for (t in seq_len(timesteps)) {
    new_grid <- grid
    
    for (i in 1:n) {
      for (j in 1:n) {
        
        state <- grid[i, j]
        neigh <- get_neighbors(grid, i, j)
        
        Pi_local <- mean(neigh == 2)
        Pj_local <- mean(neigh == 3)
        
        # ---- Occupied patches ----
        if (state == 2 && runif(1) < ei) new_grid[i, j] <- 4
        if (state == 3 && runif(1) < ej) new_grid[i, j] <- 5
        
        # ---- Latent recovery ----
        if (state == 4 && runif(1) < ri) new_grid[i, j] <- 1
        if (state == 5 && runif(1) < rj) new_grid[i, j] <- 1
        
        # ---- Colonization ----
        if (state == 1) {
          p_i <- 1 - exp(-ci * Pi_local)
          p_j <- 1 - exp(-cj * Pj_local)
          
          if (runif(1) < p_i) new_grid[i, j] <- 2
          else if (runif(1) < p_j) new_grid[i, j] <- 3
        }
        
        if (state == 4) {
          p_i <- 1 - exp(-ci * Pi_local)
          if (runif(1) < p_i) new_grid[i, j] <- 2
        }
        
        if (state == 5) {
          p_j <- 1 - exp(-cj * Pj_local)
          if (runif(1) < p_j) new_grid[i, j] <- 3
        }
      }
    }
    
    grid <- new_grid
    history[[t]] <- grid
  }
  
  history
}
time <- 250
sim <- simulate_lattice(n = 50, timesteps = time)

cols <- c("white", "blue", "red", "lightblue", "pink")
plotsteps <- round(seq(1,time,length=6))
par(mfrow = c(2, 3), mar = c(1, 1, 2, 1))
for (t in plotsteps) {
  image(
    t(sim[[t]]),
    col = cols,
    axes = FALSE,
    main = paste("t =", t)
  )
}
states <- c("empty", "Pi", "Pj", "Li", "Lj")
endstates <- mapply(function(x) {sum(sim[[t]]==x)/length(sim[[t]])}, 1:5)
names(endstates) <- states
sum(endstates)       
