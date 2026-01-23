# install.packages("deSolve")  # if needed
library(deSolve)

# ---------------------------
# Metapopulation ODE model
# ---------------------------
meta_model <- function(t, state, params) {
  with(as.list(c(state, params)), {
    
    # Recoverable empty patches
    E <- 1 - Pi - Pj - Li - Lj
    
    # ODEs
    dPi <- ci * Pi * (E + Li) - ei * Pi
    dPj <- cj * Pj * (E + Lj) - ej * Pj
    
    dLi <- ei * Pi - Li * (ri + ci * Pi)
    dLj <- ej * Pj - Lj * (rj + cj * Pj)
    
    list(c(dPi, dPj, dLi, dLj))
  })
}

# ---------------------------
# Parameters
# ---------------------------
params <- c(
  ci = 0.2,   # colonization rate of i
  cj = 0.5,   # colonization rate of j
  ei = 0.1,   # extinction rate of i
  ej = 0.1,  # extinction rate of j
  ri = 0.0001,   # recovery rate of latent i patches
  rj = 0.0001    # recovery rate of latent j patches
)

params <- c(
  ci = 0.4,   # colonization rate of i
  cj = 0.4,   # colonization rate of j
  ei = 0.1,   # extinction rate of i
  ej = 0.1,  # extinction rate of j
  ri = 0.5,   # recovery rate of latent i patches
  rj = 0.3    # recovery rate of latent j patches
)

# ---------------------------
# Initial conditions
# ---------------------------
state <- c(
  Pi = 0.1,
  Pj = 0.1,
  Li = 0.0,
  Lj = 0.0
)

# ---------------------------
# Time vector
# ---------------------------
times <- seq(0, 200, by = 0.1)

# ---------------------------
# Run simulation
# ---------------------------
out <- ode(
  y = state,
  times = times,
  func = meta_model,
  parms = params
)

out <- as.data.frame(out)

# ---------------------------
# Plot results
# ---------------------------
matplot(
  out$time,
  out[, c("Pi", "Pj", "Li", "Lj")],
  type = "l",
  lwd = 2,
  lty = 1,
  col = c("blue", "red", "blue3", "red3"),
  xlab = "Time",
  ylab = "Fraction of patches"
)

legend(
  "right",
  legend = c("P_i", "P_j", "L_i", "L_j"),
  col = c("blue", "red", "blue3", "red3"),
  lwd = 2,
  bty = "n"
)
