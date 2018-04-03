##### bound.R #####
# computes the values of the backlog bounds in the reportcd

library("dvfBm")

source("estimate_hurst.R")

# Computes the plain SNC Bound from Theorem 3.5, Equation (3.10)
# (without any statistical operations)
# n = Point in time
# x = backog
# std_dev = standard deviation,
# hurst = Hurst Parameter
# server_rate = Server Rate, also denoted C in formulas,
# arrival_rate = constant rate from the arrival model, also denoted as \lambda
backlog_bound_wrong <- function(n, x, std_dev, hurst, server_rate, arrival_rate) {
  if (server_rate <= arrival_rate) {
    stop("server rate has to be greater than the arrival rate")
  }

  k <- 1:n
  exponent <- -((x + (server_rate - arrival_rate) * k) ** 2) / (
    2 * (std_dev ** 2) * k ** (2 * hurst))
  backlog <- sum(exp(exponent))

  return(backlog)
}

# print("wrong bound:")
# print(backlog_bound(n = 10, x = 3.0, std_dev = 0.5, hurst = 0.7,
#                     server_rate = 1.0, arrival_rate = 0.6))

# Computes the plain SNC Bound from Theorem 3.10, Equation (3.12)
# (without any statistical operations)
# n = Point in time
# x = backog
# std_dev = standard deviation,
# hurst = Hurst Parameter
# server_rate = Server Rate, also denoted C in formulas,
# arrival_rate = constant rate from the arrival model, also denoted as \lambda
# tau = discretization parameter > 0
backlog_bound <- function(n, x, std_dev, hurst, server_rate,
                                arrival_rate, tau = 0.9) {
  if (server_rate <= arrival_rate) {
    stop("server rate has to be greater than the arrival rate")
  }

  k <- (floor(1 / tau) + 1):(floor(n / tau) + 1)
  exponent <- -((x - server_rate * tau + (
    server_rate - arrival_rate) * k * tau) ** 2) / (
    2 * (std_dev ** 2) * (k * tau) ** (2 * hurst))
  backlog <- sum(exp(exponent))

  return(backlog)
}

# print("discretized bound:")
# print(backlog_bound_discr(n = 10, x = 3.0, std_dev = 0.5, hurst = 0.7,
#                           server_rate = 1.0, arrival_rate = 0.6, tau = 1.0))

# for (tau in c(0.1, 0.3, 0.5, 0.7, 0.75, 0.8, 0.85, 0.9, 1.0)) {
#   print(backlog_bound_discr(n = 10, x = 3.0, std_dev = 0.5, hurst = 0.7,
#                             server_rate = 1.0, arrival_rate = 0.6, tau = tau))
# }

# numerical evuluation shows that a tau value of 0.85 is optimal for a given
# paramter set

# for (x in c(3.0, 5.0, 7.0, 10.0)) {
#   print("wrong bound:")
#   print(backlog_bound(n = 10, x = x, std_dev = 0.5, hurst = 0.7,
#                       server_rate = 1.0, arrival_rate = 0.6))
#   print("discretized bound:")
  # print(backlog_bound_discr(n = 10, x = x, std_dev = 0.5, hurst = 0.7,
  #                           server_rate = 1.0, arrival_rate = 0.6,
  #                           tau = 0.85))
  # print(backlog_bound_discr(n = 10, x = x, std_dev = 0.5, hurst = 0.7,
  #                           server_rate = 1.0, arrival_rate = 0.6,
  #                           tau = 0.9))
# }

# We observe that the discretization gap becomes less significant
# for a larger backlog
# (as the term server_rate * tau in the numerator gets smaller than
# the backlog x)


# Computes the statistical backlog bound based on the FGN increments
# (Statistical version of theorem 3.1)
# n = Point in Time
# x = Backlog
# std_dev = standard deviation
# server_rate = Server Rate, also denoted C
# arrival_rate = constant arrival rate, also denoted \lambda
# conflevel = confidence level of estimation

stat_backlog_bound <- function(flow_increments, n, x, std_dev, server_rate,
                               arrival_rate, conflevel = 0.95) {
  if (server_rate < arrival_rate) {
    stop("The server rate has to be greater than the arrival rate")
  }

  N <- length(flow_increments)

  h_estimated <- estimate_hurst(
    flow_increments = flow_increments, arrival_rate = arrival_rate,
    std_dev = std_dev)
  h_up <- conf_level_hurst(amount_increments = N, h_estimated = h_estimated,
                           conflevel = conflevel)
  # print(paste0("h_up = ", h_up))

  # TODO: find a more suitable name than "backlog_prob"
  backlog_prob <- (1 - conflevel) + backlog_bound(
    n = n, x = x, std_dev = std_dev, hurst = h_up,
    server_rate = server_rate, arrival_rate = arrival_rate)

  # print(paste0("x = ", x, ", backlog_prob = ", backlog_prob))

  return(backlog_prob)
}


# flow_example <- build_flow(arrival_rate = 1.0, hurst = 0.7, n = 2 ** 12,
#                            std_dev = 1.0)
# print(stat_backlog_bound(flow_increments = flow_example, n = 10, x = 3.0,
#                          std_dev = 1.0, server_rate = 1.0, arrival_rate = 0.6,
#                          conflevel = 0.95))

# Binary search for sufficient backlog value x s.t. P(q(n) > x) <= p,
# last parameter indicates whether SNC or stat_nc bound should be used
# n = Point in time
# p = violation probability
# std_dev = standard deviation
# hurst = Hurst parameter
# server_rate = Server Rate, also known as C in formulas
# arrival_rate = constant rate from the arrival model, also denoted \lambda
# splits = number of iterations for binary search
# conflevel = confidence level if estimation was used
# traffic = input arrival traffic (only necessary for statnc bound)
inverse_bound <- function(n, std_dev, hurst,
                          arrival_rate, server_rate, p = 10  **  (-3),
                          splits = 10, conflevel = 0.95, traffic = NaN,
                          estimate_traffic = FALSE) {
  if (estimate_traffic && is.nan(estimate_traffic)) {
    stop("We need some traffic to estimate")
  }

  if (estimate_traffic && p < (1 - conflevel)) {
    print(paste0("p = ", p, " 1 - conflevel = ", 1 - conflevel))
    stop("The bound runs in an infinite loop as the stat_backlog_bound()
         bound can never be below (1-alpha)")
  }

  backlog <- 0.5
  difference <- 1
  its <- 0

  # Search for the backlog value where bound <= p holds for the first time,
  # bisect from there

  if (estimate_traffic) {
    probbound <- stat_backlog_bound(
      flow_increments = traffic, n = n, x = backlog, std_dev = std_dev,
      server_rate = server_rate, arrival_rate = arrival_rate,
      conflevel = conflevel)
  } else {
    probbound <- backlog_bound(
      n = n, x = backlog, std_dev = std_dev, hurst = hurst,
      server_rate = server_rate, arrival_rate = arrival_rate)
  }
  while (probbound > p) {
    difference <- backlog
    backlog <- 2 * backlog
    if (estimate_traffic) {
      probbound <- stat_backlog_bound(
        flow_increments = traffic, n = n, x = backlog, std_dev = std_dev,
        server_rate = server_rate, arrival_rate = arrival_rate,
        conflevel = conflevel)
    } else {
      probbound <- backlog_bound(n = n, x = backlog, std_dev = std_dev,
                                 hurst = hurst, server_rate = server_rate,
                                 arrival_rate = arrival_rate)
    }
  }

  difference <- difference / 2
  backlog <- backlog - difference
  # Bisect $splits times
  while (its < splits) {
    if (estimate_traffic) {
      probbound <- stat_backlog_bound(
        flow_increments = traffic, n = n, x = backlog, std_dev = std_dev,
        server_rate = server_rate, arrival_rate = arrival_rate,
        conflevel = conflevel)
    } else {
      probbound <- backlog_bound(
        n = n, x = backlog, std_dev = std_dev, hurst = hurst,
        server_rate = server_rate, arrival_rate = arrival_rate)
    }
  # If the bound is smaller -> continue with "left" half, else "right"
    if (probbound <= p ) {
      difference <- difference / 2
      backlog <- backlog - difference
    } else {
      difference <- difference / 2
      backlog <- backlog + difference
    }
    its <- its + 1
  }
  return(max(0, backlog))
}

# print(inverse_bound(n = 10, p = 10  **  (-2), std_dev = 0.5, hurst = 0.7,
#                     arrival_rate = 0.6, server_rate = 1.0, splits = 10,
#                     conflevel = 0.95, traffic = NaN,
#                     estimate_traffic = FALSE))

# flow_example <- build_flow(arrival_rate = 1.0, hurst = 0.7, n = 2 ** 12,
#                            std_dev = 1.0)
# print(inverse_bound(n = 10, p = 10  **  (-2), std_dev = 0.5, hurst = 0.7,
#                     arrival_rate = 0.6, server_rate = 1.0, splits = 10,
#                     conflevel = 0.995, traffic = flow_example,
#                     estimate_traffic = TRUE))
