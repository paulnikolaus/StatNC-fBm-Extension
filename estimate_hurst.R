source("BeranWhittle.R")

# Estimates the Hurst parameter of the given traffic
# (we assume Kelly's traffic model)
# with periodograms
# flow_increments = flow_increments
# arrival_rate = constant rate of the flow, also denoted \lambda
# std_dev = standard deviation
estimate_hurst <- function(flow_increments, arrival_rate = 1.0, std_dev = 1.0) {
   # Extract the gaussian noise from increments
   N <- length(flow_increments)
   
   fgn_traffic <- (flow_increments - arrival_rate) / std_dev

   # fbm_traffic <- cumsum(fgn_traffic)

  log_frequency <- log(spec.pgram(fgn_traffic)$freq)
  log_frequency_short <- use_only_first_part(log_frequency, 0.1)

  log_periodogram <- log(spec.pgram(fgn_traffic)$spec)
  log_periodogram_short <- use_only_first_part(log_periodogram, 0.1)

  fitted <- lm(log_periodogram_short~log_frequency_short)
  # y_value <- fitted$coefficients[1]
  slope <- fitted$coefficients[2]
  
  print("slope")
  print(slope)
  h_estimated <- (1 - slope) / 2

  if (h_estimated >= 1 || h_estimated <= 0.5) {
    print("h_estimated")
    print(h_estimated)
    stop("h_estimated must be in (0.5, 1)")
  }
  
  return(h_estimated)
}

# Gives the confidence interval for the Hurst estimator
# using the periodogram approach in Rose
# amount_increments = length of the fbm_traffic vector
# h_estimated = estimated value of h
# conflevel = confidence level of the estimation

conf_level_hurst <- function(amount_increments, h_estimated,
                             conflevel = 0.95) {
  N <- amount_increments
  
  D <- CetaFGN(eta = c(H = h_estimated))
  V <- 2 * D  **  (-1)
  alpha <- 1 - conflevel
  interval <- qnorm(1 - alpha) * sqrt(V / N)
  H_up <- h_estimated + interval
  print("H_up = ")
  return(H_up)
}

use_only_first_part <- function(input_vector, share) {
  return(input_vector[1:(round(length(input_vector) * share))])
}
