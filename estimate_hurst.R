##### estimate_hurst.R #####
source("BeranWhittle.R")

# Estimates the Hurst parameter of the given traffic
# (we assume Kelly's traffic model)
# with periodograms
# TODO: Rose: In most cases, this will lead to a wrong estimate of H 
# -> better to use other estimator

# flow_increments = flow_increments
# arrival_rate = constant rate of the flow, also denoted \lambda
# std_dev = standard deviation
estimate_hurst <- function(flow_increments, arrival_rate, std_dev = 1.0) {
   # Extract the gaussian noise from flow increments
   fgn_traffic <- (flow_increments - arrival_rate) / std_dev

   # fbm_traffic <- cumsum(fgn_traffic)
   # Rose: "For FGN and FARIMA...", i.e., we don't use FBM

  log_frequency <- log(spec.pgram(fgn_traffic)$freq)
  log_frequency_short <- use_only_first_part(log_frequency, 0.1)

  log_periodogram <- log(spec.pgram(fgn_traffic)$spec)
  log_periodogram_short <- use_only_first_part(log_periodogram, 0.1)

  fitted <- lm(log_periodogram_short~log_frequency_short)
  # y_value <- fitted$coefficients[1]
  slope <- fitted$coefficients[2]
  
  # print(slope)
  h_estimated <- (1 - slope) / 2

  if (h_estimated >= 1 || h_estimated <= 0.5) {
    print("h_estimated")
    print(h_estimated)
    stop("h_estimated must be in (0.5, 1)")
  }
  
  return(h_estimated)
}

# flow_example <- build_flow(arrival_rate = 1.0, hurst = 0.7, n = 2 ** 12,
#                            std_dev = 1.0)
# print(estimate_hurst(flow_increments = flow_example, arrival_rate = 1.0,
#                      std_dev = 1.0))


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

  return(H_up)
}

# flow_example <- build_flow(arrival_rate = 1.0, hurst = 0.7, n = 2 ** 12,
#                            std_dev = 1.0)
# N <- length(flow_example)
# print(conf_level_hurst(amount_increments = N, h_estimated = 0.7,
#                        conflevel = 0.95))


use_only_first_part <- function(input_vector, share) {
  return(input_vector[1:(round(length(input_vector) * share))])
}
