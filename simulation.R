##### Simulation.R #####
## This file holds the necessary functions for the simulation.

# Libraries:
# circFbm
library("dvfBm")
# Best plots EU
library("ggplot2")

source("BeranWhittle.R") #' CetaFGN()

#' Build a fBm flow following the Kelly model.
#' First a cumulative FBM flow is built and then divided into FGN increments
#' that are easier to use in the SNC setting
#' @param rate constant arrival rate.
#' @param sample_length point in time until which the traffic
#'                      should be generated.
#' @param std_dev standard_deviation.
#' @return sample path of flow increments.
build_flow <- function(arrival_rate, hurst, sample_length, std_dev = 1.0) {
  # Previously:
  #' cumuflow <- rep(NA, sample_length)
  #' flow <- rep(NA, sample_length)
  #' fbm <- circFBM(n, hurst, plotfBm = FALSE)
  #' for (t in 1:n) {
  #'   cumuflow[t] <- arrival_rate * t + std_dev * fbm[t]
  #'   if (t >= 2) {
  #'     flow[t] <- cumuflow[t] - cumuflow[t - 1]
  #'   } else {
  #'     flow[t] <- cumuflow[t]
  #'   }
  #' }
  #' return(flow)

  # changed fbm to fbm * (sample_length ** hurst) as the package documentations
  # says that fbm is otherwise only created in the interval
  # 0, ..., (n - 1) / n
  fbm <- circFBM(n = sample_length, H = hurst, plotfBm = FALSE) * (
    sample_length ** hurst)

  # work_around to avoid completely negative traffic
  # while (mean(fbm) < 0) {
  #   fbm <- circFBM(n, hurst, plotfBm = FALSE)
  # }

  cumuflow <- arrival_rate * (1:sample_length) + std_dev * fbm
  cumuflow_shift <- c(0, cumuflow[-length(cumuflow)])

  flow_increments <- cumuflow - cumuflow_shift

  return(flow_increments)
}

#' @example
# print(build_flow(arrival_rate = 1.0, hurst = 0.7, n = 2 ** 12,
#                  std_dev = 1.0))


#' Computes the empricial backlog distribution for a specific point in time
#' for FGN arrivals with mean arrival rate, Hurst parameter hurst and
#' standard deviation std_dev at a server with the given server rate
compute_distribution <- function(
  arrival_rate, hurst, sample_length, time_n, server_rate, std_dev = 1.0,
  iterations = 10 ** 3) {
  backlogs <- rep(NA, iterations)

  for (i in 1:iterations) {
    flow <- build_flow(
      arrival_rate = arrival_rate, hurst = hurst,
      sample_length = sample_length, std_dev = std_dev)
    b <- 0
    for (k in 1:time_n) {
      # Lindley's equation:
      b <- max(b + flow[k] - server_rate, 0)
    }
    backlogs[i] <- b

    .show_progress(i, iterations, prog_msg = "compute_distribution()")
  }
  return(backlogs)
}

#' @example
# print(compute_distribution(
#   arrival_rate = 1.0, hurst = 0.7, sample_length = 10 ** 4,
#   time_n = 10 ** 4, server_rate = 2.0, std_dev = 1.0,
#   iterations = 10 ** 3))


.show_progress <- function(it, max_iterations, prog_msg = "progess") {
  perc <- round(max_iterations / 10)
  if (it %% perc == 0) {
    print(paste0(prog_msg, ": ", it * 100 / max_iterations, "%"))
  }
}
