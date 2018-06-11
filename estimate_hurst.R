##### estimate_hurst.R #####
library("fArma")
library("FGN")
# progress bar apply
library("pbapply")

source("BeranWhittle.R")
source("simulation.R")

#'Estimates the Hurst parameter of the given traffic
#'(we assume Kelly's traffic model) with periodograms
#' Rose: In most cases, this (i.e., periodogram based estimator)
#' will lead to a wrong estimate of H
#' -> better to use other estimator
# 
# Helper function for periodogram based estimation
# use_only_first_part <- function(input_vector, share) {
#   return(input_vector[1:(round(length(input_vector) * share))])
# }

#' @param flow_increments flow increments.
#' @param arrival_rate constant rate of the flow, also denoted lambda.
#' @param std_dev standard deviation, also denoted as sigma.
#' @return estimated hurst parameter.
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

#' @examples
#' flow_example <- build_flow(arrival_rate = 1.0, hurst = 0.7,
#'                            sample_length = 2 ** 14, std_dev = 1.0)
#' print(estimate_hurst(flow_increments = flow_example, arrival_rate = 1.0,
#'                      std_dev = 1.0))


#' Gives the confidence interval for the Hurst estimator
#' using the periodogram approach in Rose
#'
#' @param amount_increments length of the fbm_traffic vector.
#' @param h_estimated estimated value of h.
#' @param conflevel confidence level of the estimation.
#' @return lower and upper confidence interval as a vector.
get_h_up <- function(sample_length, h_estimated, conflevel) {
  n <- sample_length

  SD <- CetaFGN(eta = h_estimated)
  # SD <- matrix(SD, ncol = 1, nrow = 1, byrow = T) / n
  SD <- SD / n

  alpha <- (1 - conflevel)
  # we use the one-sided confidence intervall as we are only worried about
  # underestimation
  # otherwise we have to use 1 - alpha / 2
  # h_up <- h_estimated + qnorm(1 - alpha) * sqrt(SD[1, 1])
  h_up <- h_estimated + qnorm(1 - alpha) * sqrt(SD)

  # confidence interval of hurst must be in (0, 1)
  # Estimated hurst parameter cannot be above 1 as this gives an error
  if (h_up > 1.0) {
    h_up <- 1.0
  }

  return(h_up)
}

#' @example
#' flow_example <- build_flow(arrival_rate = 1.0, hurst = 0.7,
#'                            sample_length = 2 ** 12, std_dev = 1.0)
#' sample_length <- length(flow_example)
#' print(get_h_up(sample_length = sample_length, h_estimated = 0.7,
#'                conflevel = 0.95))


#' Convenience function for estimation of h's confidence interval
#'
#' @param flow_increments output of build_flow().
#' @param arrival_rate = arrival_rate used in the traffic model.
#' @param std_dev = std_dev of flow.
#' @return lower CI, estimated value, and upper CI as a vector.
flow_to_h_est_up <- function(flow_increments, arrival_rate, std_dev,
                                 conflevel) {
  sample_length <- length(flow_increments)

  h_estimated <- estimate_hurst(flow_increments = flow_increments,
                                arrival_rate = arrival_rate, std_dev = std_dev)
  h_up <- get_h_up(sample_length = sample_length, h_estimated = h_estimated,
                   conflevel = conflevel)

  return(list("h_est" = h_estimated, "h_up" = h_up))
}

flow_to_h_est_up_get_fit <- function(flow_increments, arrival_rate, std_dev) {
  # use of FGN-package

  fgn_traffic <- (flow_increments - arrival_rate) / std_dev

  res <- GetFitFGN(z = fgn_traffic, ciQ = TRUE)

  # we can only return an interval for the confidence level = 95%
  return(list("h_est" = res$"H", "h_up" = res$"ci"[2]))
}

flow_to_h_est_up_fast <- function(flow_increments, arrival_rate, std_dev) {
  # Kettani, Houssain, and John A. Gubner.
  # "A novel approach to the estimation of the Hurst parameter in
  # self-similar traffic." Local Computer Networks, 2002.
  # Proceedings. LCN 2002. 27th Annual IEEE Conference on. IEEE, 2002.
  sample_length <- length(flow_increments)
  fgn_traffic <- (flow_increments - arrival_rate) / std_dev

  rho_hat_vec <- acf(fgn_traffic, type = "correlation", plot = TRUE)
  rho_hat <- rho_hat_vec$"acf"[2]
  hurst_hat <- 0.5 * (1 + log2(1 + rho_hat))

  # we can only return an interval for the confidence level = 95%
  return(list("h_est" = hurst_hat,
              "h_up" = hurst_hat + 2.5 / sqrt(sample_length)))
}

#' @examples
# flow_example <- build_flow(arrival_rate = 1.0, hurst = 0.7,
#                            sample_length = 2 ** 11, std_dev = 1.0)
# print(flow_to_h_est_up(flow_increments = flow_example, arrival_rate = 1.0,
#                        std_dev = 1.0, conflevel = 0.95))
# print(flow_to_h_est_up_get_fit(flow_increments = flow_example,
#                                arrival_rate = 1.0, std_dev = 1.0))
# print(flow_to_h_est_up_fast(flow_increments = flow_example,
#                             arrival_rate = 1.0, std_dev = 1.0,
#                             conflevel = 0.95))


# Helper function to calculate confidence intervals
# of upper confidence interval
# TODO: delete this function
ci_help <- function(data, conf.level = 0.95) {
  # Check if all data entries are equal -> No confidence interval
  if (all(data == data[1])) {
    return(c(data[1], data[1]))
  }

  t <- t.test(data, conf.level = conf.level)$conf.int
  return(c(t[1], t[2]))
}


# Compute a confidence interval for the estimation of H
# TODO: delete this function
# confint_of_h_up <- function(
#   sample_length, arrival_rate, hurst, std_dev, conflevel, iterations,
#   confint.conflevel) {
#   hurst_up_estimates <- rep(NA, iterations)
#   for (i in 1:iterations) {
#     f <- build_flow(
#       arrival_rate = arrival_rate, hurst = hurst,
#       sample_length = sample_length, std_dev = std_dev)
#     hurst_up_estimates[i] <- flow_to_h_est_up(
#       flow_increments = f, arrival_rate = arrival_rate,
#       std_dev = std_dev, conflevel = conflevel)$"h_up"
#     # hurst_up_estimates[i] <- flow_to_h_est_up(
#     #   flow_increments = f, arrival_rate = arrival_rate, std_dev = std_dev,
#     #   conflevel = conflevel)$"h_up"
# 
#     .show_progress(i, iterations, prog_msg = "confint_of_h_up()")
#   }
#   ci <- ci_help(data = hurst_up_estimates, conf.level = confint.conflevel)
#   m <- mean(hurst_up_estimates)
# 
#   return(append(m, ci))
# }

# print(confint_of_h_up(
#   sample_length = 2 ** 12, arrival_rate = 1.0, hurst = 0.7, std_dev = 1.0,
#   conflevel = 0.999, iterations = 10 ** 2, confint.conflevel = 0.999))


#' Compute mean of the confidence interval's upper value
#' @return vector of estimated h_up's.
est_h_up_vector <- function(
  sample_length, arrival_rate, hurst, std_dev, conflevel, iterations) {
  # old version with for-loop:
  # hurst_up_estimates <- rep(NA, iterations)
  # for (i in 1:iterations) {
  #   f <- build_flow(
  #     arrival_rate = arrival_rate, hurst = hurst,
  #     sample_length = sample_length, std_dev = std_dev)
  #   hurst_up_estimates[i] <- flow_to_h_est_up(
  #     flow_increments = f, arrival_rate = arrival_rate, std_dev = std_dev,
  #     conflevel = conflevel)$"h_up"
  # }

  # added input parameter in order to use sapply()
  build_flow_iter <- function(iter) {
    return(build_flow(
      arrival_rate = arrival_rate, hurst = hurst,
      sample_length = sample_length, std_dev = std_dev))
  }

  flow_to_h_up <- function(flow_increments) {
    return(flow_to_h_est_up(
      flow_increments = flow_increments, arrival_rate = arrival_rate,
      std_dev = std_dev, conflevel = conflevel)$"h_up")
  }

  # flow_to_h_est <- function(flow_increments) {
  #   return(flow_to_h_est_up(
  #     flow_increments = flow_increments, arrival_rate = arrival_rate,
  #     std_dev = std_dev, conflevel = conflevel)$"h_est")
  # }

  flow_matrix <- sapply(1:iterations, build_flow_iter)
  # dim(flowmatrix) = sample_length  iterations

  #' @example
  #' print("hurst_estimates")
  #' print(mean(pbapply(flow_matrix, 2, flow_to_h_est)))
  #' result: 0.7032785

  #' print("hurst_up_estimates")
  #' print(mean(pbapply(flow_matrix, 2, flow_to_h_up)))
  #' result: 0.7112128

  hurst_up_estimates <- pbapply(flow_matrix, 2, flow_to_h_up)
  # pbapply() = apply() with progress bar

  return(hurst_up_estimates)
}


#' Helper function. Takes a vector of repeatedly estimated Hurst parameters
#' and returns the lower and upper quantile.
#' @return lower quantile of h_ups, mean of h_up, upper quantile of h_up
compute_h_up_quantile <- function(h_vector, quantile_prob = 0.95) {
  hurst_up_means <- mean(h_vector)

  beta <- 1 - quantile_prob
  return(list("Hurst_lower_quant" = quantile(h_vector, beta / 2)[[1]],
              "Hurst_up_mean" = hurst_up_means,
              "Hurst_upper_quant" = quantile(h_vector, 1 - beta / 2)[[1]]))

  # h_confidenceInterval = ci_help(h_vector)
  # return(list("Hurst_lower_quant" = h_confidenceInterval[[1]],
  #             "Hurst_up_mean" = hurst_up_means,
  #             "Hurst_upper_quant" = h_confidenceInterval[[2]]))

}

#' @examples 
# h_ups <- est_h_up_vector(sample_length = 2 ** 13, arrival_rate = 1.0,
#                          hurst = 0.7, std_dev = 1.0, conflevel = 0.999,
#                          iterations = 100)
# print(compute_h_up_quantile(h_vector = h_ups))
