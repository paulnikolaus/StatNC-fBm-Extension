##### estimate_hurst.R #####
library("fArma")

source("BeranWhittle.R")
source("simulation.R")

# Estimates the Hurst parameter of the given traffic
# (we assume Kelly's traffic model)
# with periodograms
# Rose: In most cases, this (i.e., periodogram based estimator)
# will lead to a wrong estimate of H
# -> better to use other estimator

# Helper function for periodogram based estimation
# use_only_first_part <- function(input_vector, share) {
#   return(input_vector[1:(round(length(input_vector) * share))])
# }

# flow_increments = flow_increments
# arrival_rate = constant rate of the flow, also denoted \lambda
# std_dev = standard deviation
estimate_hurst <- function(flow_increments, arrival_rate, std_dev = 1.0) {
   # Extract the gaussian noise from flow increments
   fgn_traffic <- (flow_increments - arrival_rate) / std_dev

  # old, self-written, periodogram approach

  # log_frequency <- log(spec.pgram(fgn_traffic, plot = FALSE)$freq)
  # log_frequency_short <- use_only_first_part(log_frequency, 0.1)
  #
  # log_periodogram <- log(spec.pgram(fgn_traffic, plot = FALSE)$spec)
  # log_periodogram_short <- use_only_first_part(log_periodogram, 0.1)
  #
  # fitted <- lm(log_periodogram_short~log_frequency_short)
  # # y_value <- fitted$coefficients[1]
  # slope <- fitted$coefficients[2]
  #
  # # print(slope)
  # h_estimated <- (1 - slope) / 2


  h_estimated <- perFit(x = fgn_traffic)@hurst$"H"

  # use str(class) to see all possibilities to obtain values of
  #  an S4 class

  if (h_estimated >= 1 || h_estimated <= 0.5) {
    # print("h_estimated")
    # print(h_estimated)
    warning("h_estimated must be in (0.5, 1)")
  }

  return(h_estimated)
}

flow_example <- build_flow(arrival_rate = 1.0, hurst = 0.7,
                           sample_length = 2 ** 14, std_dev = 1.0)
print(estimate_hurst(flow_increments = flow_example, arrival_rate = 1.0,
                     std_dev = 1.0))


# Gives the confidence interval for the Hurst estimator
# using the periodogram approach in Rose
# amount_increments = length of the fbm_traffic vector
# h_estimated = estimated value of h
# conflevel = confidence level of the estimation

get_h_up <- function(amount_increments, h_estimated,
                             conflevel = 0.95) {
  N <- amount_increments

  D <- CetaFGN(eta = c(H = h_estimated))
  V <- 2 * D  **  (-1)
  alpha <- 1 - conflevel
  h_up <- h_estimated + qnorm(1 - alpha) * sqrt(V / N)

  # confidence interval of hurst must be in (0, 1)
  # Estimated hurst parameter cannot be above 1 as this gives an error
  if (h_up > 1.0) {
    h_up <- 1.0
  }

  return(h_up)
}

# flow_example <- build_flow(arrival_rate = 1.0, hurst = 0.7,
#                            sample_length = 2 ** 12, std_dev = 1.0)
# amount_increments <- length(flow_example)
# print(get_h_up(amount_increments = amount_increments,
#                        h_estimated = 0.7, conflevel = 0.95))


# Convenience function for estimation of h_up
# flow_increments = FGN Flow
# arrival_rate = arrival_rate used in the traffic model
# std_dev = std_dev of flow

flow_to_h_up <- function(flow_increments, arrival_rate, std_dev, conflevel) {
  amount_increments <- length(flow_increments)

  h_estimated <- estimate_hurst(
    flow_increments = flow_increments, arrival_rate = arrival_rate,
    std_dev = std_dev)
  h_up <- get_h_up(amount_increments = amount_increments,
                   h_estimated = h_estimated, conflevel = conflevel)
  # print(paste0("h_up = ", h_up))
  return(h_up)
}


# Helper function to calculate confidence intervals
# of upper confidence interval
ci_help <- function(data, conf.level = 0.95) {
  # Check if all data entries are equal -> No confidence interval
  if (all(data == data[1])) {
    return(c(data[1], data[1]))
  }

  t <- t.test(data, conf.level = conf.level)$conf.int
  return(c(t[1], t[2]))
}


# Compute a confidence interval for the estimation of H
# conflevel: confidence level of hurst estimation
# confint.conflevel: confidence level for upper hurst confidence level
confint_of_h_up <- function(
  sample_length, arrival_rate, hurst, std_dev, conflevel, iterations,
  confint.conflevel) {
  hurst_up_estimates <- rep(NA, iterations)
  for (i in 1:iterations) {
    f <- build_flow(
      arrival_rate = arrival_rate, hurst = hurst,
      sample_length = sample_length, std_dev = std_dev)
    hurst_up_estimates[i] <- flow_to_h_up(
      flow_increments = f, arrival_rate = arrival_rate, std_dev = std_dev,
      conflevel = conflevel)
  }
  ci <- ci_help(data = hurst_up_estimates, conf.level = confint.conflevel)
  m <- mean(hurst_up_estimates)
  return(append(m, ci))
}
# print(confint_of_h_up(
#   sample_length = 2 ** 12, arrival_rate = 1.0, hurst = 0.7, std_dev = 1.0,
#   conflevel = 0.999, iterations = 10 ** 2, confint.conflevel = 0.999))


# Compute mean of the confidence interval's upper value

mean_of_h_up <- function(
  sample_length, arrival_rate, hurst, std_dev, conflevel, iterations) {
  # old version with for-loop:
  # hurst_up_estimates <- rep(NA, iterations)
  # for (i in 1:iterations) {
  #   f <- build_flow(
  #     arrival_rate = arrival_rate, hurst = hurst,
  #     sample_length = sample_length, std_dev = std_dev)
  #   hurst_up_estimates[i] <- flow_to_h_up(
  #     flow_increments = f, arrival_rate = arrival_rate, std_dev = std_dev,
  #     conflevel = conflevel)
  # }

  build_flow_iter <- function(
    iter, arrival_rate = arrival_rate, hurst = hurst,
    sample_length = sample_length, std_dev = std_dev) {
    return(build_flow(
      arrival_rate = arrival_rate, hurst = hurst,
      sample_length = sample_length, std_dev = std_dev))
  }

  flow_matrix <- sapply(1:iterations, build_flow_iter,
                        arrival_rate = arrival_rate, hurst = hurst,
                        sample_length = sample_length, std_dev = std_dev)
  # dim(flowmatrix) = sample_length  iterations
  hurst_up_estimates <- apply(flow_matrix, 2, flow_to_h_up,
                              arrival_rate = arrival_rate, std_dev = std_dev,
                              conflevel = conflevel)


  return(mean(hurst_up_estimates))
}

# print(mean_of_h_up(
#   sample_length = 2 ** 12, arrival_rate = 1.0, hurst = 0.7, std_dev = 1.0,
#   conflevel = 0.999, iterations = 10 ** 2))
