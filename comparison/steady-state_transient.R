##### steady-state_transient.R #####

# Computes the plain (not discretized) SNC Bound from Theorem 3.5,
# Equation (3.10)
# (without any statistical operations)
# time_n = Point in time
# x = backog
# std_dev = standard deviation,
# hurst = Hurst Parameter
# server_rate = Server Rate, also denoted C in formulas,
# arrival_rate = constant rate from the arrival model, also denoted as \lambda
backlog_bound_transient <- function(time_n, x, std_dev, hurst, server_rate,
                                    arrival_rate) {
  if (server_rate <= arrival_rate) {
    stop("server rate has to be greater than the arrival rate")
  }

  k <- 1:time_n
  exponent <- -((x + (server_rate - arrival_rate) * k) ** 2) / (
    2 * (std_dev ** 2) * k ** (2 * hurst))
  backlog <- sum(exp(exponent))

  return(backlog)
}

backlog_bound_steady <- function(std_dev, hurst, server_rate, arrival_rate) {
  if (server_rate <= arrival_rate) {
    stop("server rate has to be greater than the arrival rate")
  }
  numerator <- gamma(1 / (2 * (1 - hurst)))
  denominator <- 2 * (1 - hurst) * (
    0.5 * ((server_rate - arrival_rate) ** 2) / (std_dev ** 2)) ** (
      0.5 / (1 - hurst))

  return(numerator / denominator)
}

sigma <- 1.0
h <- 0.7
C <- 3.0
lamb <- 1.0

print("transient bound:")
print(backlog_bound_transient(time_n = 10, x = 3.0, std_dev = sigma, hurst = h,
                              server_rate = C, arrival_rate = lamb))
print("steady-state bound:")
print(backlog_bound_steady(std_dev = sigma, hurst = h, server_rate = C,
                           arrival_rate = lamb))

#"transient bound:"
#0.0008727755
#"steady-state bound:"
#0.4739116

# We observe that the steady-state bound is about a 1000 times worse
