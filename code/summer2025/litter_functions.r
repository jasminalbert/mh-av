#update_litter <- function(N1, N2, L, litter_rate1 = 0.01, litter_rate2 = 0.01, decay = 0.05) {
  # Inputs:
  #   N1, N2: species biomass matrices
  #   L: current litter matrix
  #   litter_rate1, litter_rate2: rate at which each species produces litter
  #   decay: fraction of litter that decays each time step

#  litter_input <- litter_rate1 * N1 + litter_rate2 * N2
#  L_new <- (1 - decay) * L + litter_input

#  return(L_new)
#}


# Convert (i,j) to 1D index
ij_to_index <- function(i, j, nrows) {
  (j - 1) * nrows + i
}

# Initialize: flat list of cohort lists
initialize_litter_cohorts <- function(grid_size) {
  litter_cohorts <- vector("list", grid_size^2)
  lapply(litter_cohorts, function(X){data.frame(matrix(NA,nrow=0,ncol=3, dimnames=list(NULL,c("amount1", "amount2","birth_time"))))})
}

# Add litter cohort
add_litter <- function(litter_cohorts, N1, N2, time_step) {
  new_litter1 <- 0.1*N1
  #print(new_litter1)
  new_litter2 <- 0.2*N2
  #print(new_litter2)
  new_litter <- new_litter1+new_litter2
  #print(new_litter)
  n <- nrow(new_litter1)
  for (i in 1:n) {
    for (j in 1:n) {
      idx <- ij_to_index(i, j, n)
      if (new_litter[i, j] > 0) {
        litter_cohorts[[idx]] <- rbind(
          litter_cohorts[[idx]],
          data.frame(amount1 = new_litter1[i, j], amount2 = new_litter2[i, j], birth_time = time_step))
      }
    }
  }
  litter_cohorts
}

# Decay litter and return total matrix + updated cohorts
decay_litter <- function(litter_cohorts, current_time, decay_rate1, decay_rate2, grid_size) {
  total_litter <- matrix(0, grid_size, grid_size)
  for (i in 1:grid_size) {
    for (j in 1:grid_size) {
      idx <- ij_to_index(i, j, grid_size)
      cohorts <- litter_cohorts[[idx]]
      updated_cohorts <- matrix(nrow=0,ncol=2)
      for (c in 1:nrow(cohorts)) {
        age <- current_time - cohorts[c,]$birth_time
        remaining1 <- cohorts[c,]$amount1 * exp(-decay_rate1 * age)
        remaining2 <- cohorts[c,]$amount2 * exp(-decay_rate2 * age)
        if (sum(remaining1,remaining2) > 1e-6) {
          updated_cohorts <- rbind(updated_cohorts, cohorts[c,])  # Keep original birth_time
        }
        total_litter[i, j] <- total_litter[i, j] + remaining1+remaining2
      }
      litter_cohorts[[idx]] <- updated_cohorts
    }
  }
  list(litter = total_litter, updated_cohorts = litter_cohorts)
}




