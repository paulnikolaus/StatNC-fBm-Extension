##### Simulation.R #####
## This file holds the necessary functions for the simulation.

# Libraries.
# For circFbm
library("dvfBm")
# Best plots EU
library("ggplot2")

source("BeranWhittle.R") # need CetaFGN()
source("Bound.R") # inverse_bound()

# Build a fBm flow following the Kelly model.
# First a cumulative FBM flow is built and then divided into FGN increments
# which are easier to use in the SNC setting
# rate = constant arrival rate, n= point in time until
# which the traffic shall be generated, std_dev = standard_deviation.

build_flow <- function(arrival_rate = 1.0, H = 0.7, n = 2^12, std_dev = 1.0) {
  fbm <- circFBM(n, H, plotfBm = FALSE)

  cumuflow <- arrival_rate * (1:n) + std_dev * fbm
  cumuflow_shift <- c(0, cumuflow[-length(cumuflow)])

  flow <- cumuflow - cumuflow_shift

  return(flow)
}

# Given a flow and a constant rate for the server, compute the backlog at
# each point in time
# Flow = FGN Arrival Flow, server_rate = constant Server rate, n=point in time
# until which the backlog shall be simulated

simulate_system <- function(flow, server_rate = 2.0, n = 2^12) {
  backlog <- rep(NA, n)
  backlog[1] <- 0
  for (i in 2:n) {
  	backlog[i] <- max(backlog[i - 1] + flow[i - 1] - server_rate, 0)
  }
  return(backlog)
}

# Computes the empricial backlog distribution for a specific point in time
# for FGN arrivals with mean arrival rate, Hurst parameter H and
# standard deviation std_dev at a server with the given server rate

compute_distribution <- function(iterations = 10^6, arrival_rate = 1.0,
                                 H = 0.7, n = 10^4, std_dev = 1.0,
                                 server_rate = 2.0) {
  backlogs <- rep(NA, iterations)

  for (i in 1:iterations) {
    flow <- build_flow(arrival_rate = arrival_rate, H = H, n = n,
                       std_dev = std_dev)
    b <- 0
    for (k in 1:n) {
      b <- max(b + flow[k] - server_rate, 0)
    }
    .showProgress(i, iterations)
    backlogs[i] <- b
  }
  return(backlogs)
}

.showProgress <- function(it, max_iterations) {
  perc <- round(max_iterations / 10)
  if (it %% perc == 0) {
    print(it * 100 / max_iterations)
  }
}
