##### bound.R #####
# computes the values of the backlog bounds in the reportcd

library("dvfBm")

source("estimate_hurst.R")
source("simulation.R")

#' Computes the plain SNC Bound from Theorem 3.10, Equation (3.12)
#' (without any statistical operations)
#' @param time_n Point in time.
#' @param x backlog.
#' @param std_dev standard deviation.
#' @param hurst Hurst Parameter.
#' @param server_rate Server Rate, also denoted C in formulas.
#' @param arrival_rate constant rate from the arrival model, also
#' denoted as lambda
#' @param tau discretization parameter > 0.
#' @return backlog violation probability.
backlog_bound <- function(time_n, x, std_dev, hurst, server_rate,
                          arrival_rate, tau = 0.9) {
  # TODO: Better default value for tau

  if (server_rate <= arrival_rate) {
    stop("server rate has to be greater than the arrival rate")
  }
  if (x - arrival_rate * tau <= 0) {
    warning("theta's sign constraint is violated")
  }

  k <- 1:(floor(time_n / tau) + 1)

  exponent <- -((x - server_rate * tau + (
    server_rate - arrival_rate) * k * tau)**2) / (
    2 * (std_dev**2) * (k * tau)**(2 * hurst))

  backlog <- sum(exp(exponent))

  return(backlog)
}

#' @examples
#' print("backlog bound:")
#' print(backlog_bound(time_n = 10, x = 3.0, std_dev = 0.5, hurst = 0.7,
#'                     server_rate = 1.0, arrival_rate = 0.6, tau = 1.0))
#
#' for (tau in c(0.1, 0.3, 0.5, 0.7, 0.75, 0.8, 0.85, 0.9, 1.0)) {
#'   print(paste0("tau: ", tau))
#'   print(paste0("bound: ", backlog_bound(
#'     time_n = 10, x = 3.0, std_dev = 0.5, hurst = 0.7,
#'     server_rate = 1.0, arrival_rate = 0.6, tau = tau)))
#' }

#' numerical evuluation shows that a tau value of 0.85 is optimal for this
#' paramter set

#' for (x in c(3.0, 5.0, 7.0, 10.0)) {
#'   print("backlog bound:")
#' print(backlog_bound(time_n = 10, x = x, std_dev = 0.5, hurst = 0.7,
#'                           server_rate = 1.0, arrival_rate = 0.6,
#'                           tau = 0.85))
#' print(backlog_bound(time_n = 10, x = x, std_dev = 0.5, hurst = 0.7,
#'                           server_rate = 1.0, arrival_rate = 0.6,
#'                           tau = 0.9))
#' }


#' Computes the statistical backlog bound based on the FGN increments
#' (Statistical version of theorem 3.1)
#' @param time_n Point in Time.
#' @param x Backlog.
#' @param std_dev standard deviation.
#' @param hurst the estimated hurst parameter.
#' @param server_rate Server Rate, also denoted C.
#' @param arrival_rate constant arrival rate, also denoted as lambda.
#' @param conflevel confidence level of estimation.
#' @return StatNC backlog violation probability.
stat_backlog_bound <- function(time_n, x, std_dev, hurst, server_rate,
                               arrival_rate, conflevel = 0.95) {
  if (server_rate < arrival_rate) {
    stop("The server rate has to be greater than the arrival rate")
  }

  backlog_stat <- (1 - conflevel) + backlog_bound(
    time_n = time_n, x = x, std_dev = std_dev, hurst = hurst,
    server_rate = server_rate, arrival_rate = arrival_rate
  )

  return(backlog_stat)
}

#' @example
#' print(stat_backlog_bound(time_n = 100, x = 3.0, std_dev = 1.0, hurst = 0.7,
#'                          server_rate = 1.0, arrival_rate = 0.6,
#'                          conflevel = 0.95))


#' Binary search for sufficient backlog value x s.t. P(q(n) > x) <= p,
#' last parameter indicates whether SNC or stat_nc bound should be used
#' @param time_n Point in time.
#' @param p violation probability.
#' @param std_dev standard deviation.
#' @param hurst Hurst parameter.
#' @param server_rate Server Rate, also known as C in formulas.
#' @param arrival_rate constant rate from the arrival model, also
#' denoted as lambda.
#' @param splits number of iterations for binary search.
#' @param conflevel confidence level if estimation was used.
#' @return backlog bound for a given probability.
inverse_bound <- function(time_n, std_dev, hurst,
                          arrival_rate, server_rate, p = 10**(-3),
                          splits = 10, conflevel = 0.95,
                          estimated_h = FALSE) {
  if (estimated_h && p < (1 - conflevel)) {
    stop(paste0("p = ", p, " < (1 - conflevel) = ", 1 - conflevel, ". \n
The bound runs in an infinite loop as the stat_backlog_bound() bound can never
                be below (1-alpha)"))
  }

  backlog <- 0.5
  difference <- 1
  its <- 0

  stat_backlog_bound_short <- function(backlog) {
    return(stat_backlog_bound(
      time_n = time_n, x = backlog, std_dev = std_dev, hurst = hurst,
      server_rate = server_rate, arrival_rate = arrival_rate,
      conflevel = conflevel
    ))
  }

  backlog_bound_short <- function(backlog) {
    return(backlog_bound(
      time_n = time_n, x = backlog, std_dev = std_dev, hurst = hurst,
      server_rate = server_rate, arrival_rate = arrival_rate
    ))
  }

  #' Search for the backlog value where bound <= p holds for the first time,
  #' bisect from there

  if (estimated_h) {
    probbound <- stat_backlog_bound_short(backlog = backlog)
  } else {
    probbound <- backlog_bound_short(backlog = backlog)
  }
  while (probbound > p) {
    difference <- backlog
    backlog <- 2 * backlog
    if (estimated_h) {
      probbound <- stat_backlog_bound_short(backlog = backlog)
    } else {
      probbound <- backlog_bound_short(backlog = backlog)
    }
  }

  difference <- difference / 2
  backlog <- backlog - difference
  #' Bisect $splits times
  while (its < splits) {
    if (estimated_h) {
      probbound <- stat_backlog_bound_short(backlog = backlog)
    } else {
      probbound <- backlog_bound_short(backlog = backlog)
    }
    #' If the bound is smaller -> continue with "left" half, else "right"
    difference <- difference / 2
    if (probbound <= p) {
      backlog <- backlog - difference
    } else {
      backlog <- backlog + difference
    }
    its <- its + 1
  }
  return(max(0, backlog))
}

#' @examples
#' print(inverse_bound(time_n = 100, p = 10  **  (-2), std_dev = 0.5,
#'                     hurst = 0.7, arrival_rate = 0.6, server_rate = 1.0,
#'                     splits = 10, conflevel = 0.999, estimated_h = FALSE))
# # result: 13.86328
#
#' flow_example <- build_flow(arrival_rate = 1.0, hurst = 0.7,
#'                            sample_length = 2 ** 12, std_dev = 1.0)
#' print(inverse_bound(time_n = 100, p = 10  **  (-2), std_dev = 0.5,
#'                     hurst = 0.7, arrival_rate = 0.6, server_rate = 1.0,
#'                     splits = 10, conflevel = 0.999, estimated_h = TRUE))
# # result: 14.08984
