#simulation of q(N)
# Load necessary library
library(dplyr)

# Function to generate synthetic data
generate_synthetic_data <- function(N, T) {
  sort(runif(N, min = 0, max = T))
}

compute_lrs <- function(data, T0){
  N <- length(data)
  N_T0 <- sum(data <=T0)
  N_1_T0 <- N-N_T0
  if (N_T0 == 0 || N_1_T0 == 0) return(-Inf)
  log_likelihood_0 <- 2 * N * log(N)
  log_likelihood_1 <- 2 * (N_T0 * log(N_T0 / T0) + N_1_T0 * log((N - N_T0) / (1 - T0)))
  
  # Likelihood ratio statistic
  statistic <- log_likelihood_1 - log_likelihood_0
  return(statistic)
}

monte_carlo_simulation <- function(M, N, T) {
  max_likelihood_ratios <- numeric(M)
  
  for (i in 1:M) {
    data <- generate_synthetic_data(N, T)
    likelihood_ratios <- sapply(data, function(T0) compute_lrs(data, T0))
    max_likelihood_ratios[i] <- max(likelihood_ratios)
  }
  
  q_N <- (mean(max_likelihood_ratios)) - 1/M
  return(q_N)
}

M <- 10000
N <- 3457  # Number of data points
N_1 <- 2377
N_2 <- 1531
T <- 1     # Interval [0, T]
# Run Monte Carlo simulation
q_N <- monte_carlo_simulation(M, N, T)
## 6.225553
q_N_1 <- monte_carlo_simulation(M, N_1, T)
## 6.061035
q_N_2 <- monte_carlo_simulation(M, N_2, T)
## 5.86771

