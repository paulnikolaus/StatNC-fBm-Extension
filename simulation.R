##### Simulation.R #####
## This file holds the necessary functions for the simulation.

# Libraries.
# For circFbm
library("dvfBm")
# Best plots EU
library("ggplot2")

source("BeranWhittle.R") # CetaFGN()
source("Bound.R") # inverse_bound()

# Build a fBm flow following the Kelly model.
# First a cumulative FBM flow is built and then divided into FGN increments
# which are easier to use in the SNC setting
# rate = constant arrival rate, n = point in time until
# which the traffic shall be generated, std_dev = standard_deviation.

build_flow <- function(arrival_rate, hurst, n, std_dev = 1.0) {
  # Previously:
  # cumuflow <- rep(NA, n)
  # flow <- rep(NA, n)
  # fbm <- circFBM(n, hurst, plotfBm = FALSE)
  # for (t in 1:n) {
  #   cumuflow[t] <- arrival_rate * t + std_dev * fbm[t]
  #
  #   if (t >= 2) {
  #     flow[t] <- cumuflow[t] - cumuflow[t - 1]
  #   } else {
  #     flow[t] <- cumuflow[t]
  #   }
  # }
  #
  # return(flow)

  # changed fbm to fbm * (n ** hurst) as the package documentations
  # says that fbm is otherwise only created in the interval
  # 0, ..., (n - 1) / n
  fbm <- circFBM(n = n, H = hurst, plotfBm = FALSE) * (n ** hurst)

  # work_around to avoid completely negative traffic
  # while (mean(fbm) < 0) {
  #   fbm <- circFBM(n, hurst, plotfBm = FALSE)
  # }

  cumuflow <- arrival_rate * (1:n) + std_dev * fbm
  cumuflow_shift <- c(0, cumuflow[-length(cumuflow)])

  flow_increments <- cumuflow - cumuflow_shift

  return(flow_increments)
}

# print(build_flow(arrival_rate = 1.0, hurst = 0.7, n = 2 ** 12,
#                  std_dev = 1.0))


# Given a flow and a constant rate for the server, compute the backlog at
# each point in time
# Flow = FGN Arrival Flow, server_rate = constant Server rate, n=point in time
# until which the backlog should be simulated

simulate_system <- function(flow_increments, server_rate, n) {
  backlog <- rep(NA, n)
  backlog[1] <- 0
  for (i in 2:n) {
    backlog[i] <- max(backlog[i - 1] + flow_increments[i - 1] - server_rate, 0)
  }
  return(backlog)
}

# flow_example <- build_flow(arrival_rate = 1.0, hurst = 0.7, n = 2 ** 12,
#                            std_dev = 1.0)
# print(simulate_system(flow_increments = flow_example, server_rate = 2.0,
#                       n = 2 ** 12))


# Computes the empricial backlog distribution for a specific point in time
# for FGN arrivals with mean arrival rate, Hurst parameter hurst and
# standard deviation std_dev at a server with the given server rate

compute_distribution <- function(arrival_rate, hurst, n, server_rate,
                                 std_dev = 1.0, iterations = 10 ** 3) {
  backlogs <- rep(NA, iterations)

  for (i in 1:iterations) {
    flow <- build_flow(arrival_rate = arrival_rate, hurst = hurst, n = n,
                       std_dev = std_dev)
    b <- 0
    for (k in 1:n) {
      b <- max(b + flow[k] - server_rate, 0)
    }
    .show_progress(i, iterations)
    backlogs[i] <- b
  }
  return(backlogs)
}

# Helper function to calculate confidence intervals
CI = function(data, conf.level = 0.95) {
  # Check if all data entries are equal -> No confidence interval
  if(all(data == data[1])) {
    return(c(data[1], data[1]))
  }
  
  t = t.test(data, conf.level = conf.level)$conf.int
  return(c(t[1], t[2]))
}

# Compute a confidence interval for the estimation of H
confint_h_up <- function(arrival_rate, hurst, std_dev, conflevel, iterations, confint.conflevel) {
  hurst_estimates = rep(NA, iterations)
  for (i in 1:iterations) {
    f <- build_flow(
      arrival_rate = arrival_rate, hurst = hurst, n = 2 ** 10, std_dev = std_dev)
    hurst_estimates[i] <- flow_to_h_up(f, arrival_rate = arrival_rate, std_dev = std_dev, conflevel = conflevel)
  }
  ci = CI(hurst_estimates, conf.level = confint.conflevel)
  m = mean(hurst_estimates)
  return(append(m, ci))
}

# print(compute_distribution(arrival_rate = 1.0, hurst = 0.7, n = 10 ** 4,
#                            server_rate = 2.0, std_dev = 1.0,
#                            iterations = 10 ** 3))


.show_progress <- function(it, max_iterations) {
  perc <- round(max_iterations / 10)
  if (it %% perc == 0) {
    print(paste0(it * 100 / max_iterations, "%"))
  }
}
