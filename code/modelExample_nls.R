# Load required library for nonlinear fitting (using nls this time)
# library(minpack.lm) # Not needed now since we're using nls directly

# Define the Beverton-Holt function for one-step projection
beverton_holt <- function(N_t, R0, alpha) {
  return((R0 * N_t) / (1 + alpha * N_t))
}

# Example data: time points and observed population sizes
# Replace 'time' and 'observed_pop' with your actual data
time <- 1:10  # Time points
observed_pop <- c(50, 70, 85, 80, 75, 68, 60, 55, 52, 50)  # Observed population sizes

# Prepare the data for fitting
N_t <- observed_pop[-length(observed_pop)]
N_tp1 <- observed_pop[-1]  # Population at time t+1

# Fit the model using non-linear least squares (nls)
# Replace 'R0_start' and 'alpha_start' with initial parameter guesses
R0_start <- 1.5  # Initial guess for R0
alpha_start <- 0.02  # Initial guess for alpha

# Define the model fitting function
fit_model <- nls(
  N_tp1 ~ (R0 * N_t) / (1 + alpha * N_t),
  start = list(R0 = R0_start, alpha = alpha_start)
)

# Output the fitting summary
summary(fit_model)

# Extract fitted parameters
fitted_params <- coef(fit_model)
R0 <- fitted_params["R0"]
alpha <- fitted_params["alpha"]

# Print the fitted parameters
cat("Fitted Parameters:\n")
cat("R0 (Reproductive Rate):", R0, "\n")
cat("Alpha (Density Dependence):", alpha, "\n")

# Generate fitted values using the fitted parameters
fitted_values <- numeric(length(time))
fitted_values[1] <- observed_pop[1]  # Start with the initial observed population
for (i in 2:length(time)) {
  fitted_values[i] <- beverton_holt(fitted_values[i - 1], R0, alpha)
}

# Plot observed vs. fitted values
plot(time, observed_pop, pch = 16, col = 'blue', ylim = range(c(observed_pop, fitted_values)),
     xlab = "Time", ylab = "Population Size", main = "Beverton-Holt Model Fit")
lines(time, fitted_values, col = 'red', lwd = 2)
legend("topright", legend = c("Observed", "Fitted"), col = c("blue", "red"), pch = c(16, NA), lty = c(NA, 1))
