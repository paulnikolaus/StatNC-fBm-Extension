##### plot__false_SNC.R #####

# This file holds the necessary functions for plotting the "wrong"
# behavior of normal StatNC.
# To this end, we include the non-parametric estimator from the StatNC paper,
# as well as the exp. one.

library("ggplot2")

source("Bound.R") # inverse_bound(), loads estimate_hurst.R and simulation.R

#' Estimate a lower bound on the distribution parameter lambda, for
#' exp. i.i.d. arrivals
#' @param  arrivals the arrivals whose parameter is estimated
#' (increments not cumulative).
#' @param conflevel confidence level of the estimator.

exp_traffic_estimator <- function(arrivals, conflevel) {
  num_samples <- length(arrivals)
  sample_sum <- sum(arrivals)
  # (non-central) Chi-Squared Distribution
  return((qchisq(conflevel, df = 2 * num_samples)) / (2 * sample_sum))
  # = lambda_min
}

#' Non parametric estimator, returns the estimated mgf, M_ is bandwidth
#' limitation on traffic increments
#' @param arrivals the arrivals whose parameter is estimated
#' (increments not cumulative).
#' @param conflevel confidence level of the estimator.
#' @param blimit bandwidth limitation.

nonparametric_estimator <- function(arrivals, conflevel, blimit) {
  num_samples <- length(arrivals)
  sample_mean_theta <- function(t) sum(exp(t * arrivals)) / num_samples
  constant <- sqrt((-log(conflevel / 2)) / (2 * num_samples))
  return(
    function(m, n, t) (sample_mean_theta(t) + constant * (
      exp(t * blimit) - 1)) ** (n - m))
}

#' Theta-dependent inverse backlog bound. Gets an estimated mgf and
#' outputs a backlog bound.
#' Has to be optimizedin theta. NOTE: Actual Viol. Prob is p - alpha
#' @param viol_prob target violation probaility.
#' @param t_time point in time which should be evaluated.
#' @param conflevel confidence level of the estimator (alpha in the paper).
#' @param emgf the estimated mgf.
#' @param server_rate constant rate of the server.
#' @param theta free variable of every mgf, has to be optimized later.
statNC_inverse_backlog_theta <- function(
  viol_prob, t_time, conflevel, emgf, server_rate, theta) {
  if (viol_prob <= 1 - conflevel) {
    stop("confidence level has to be smaller than violation probability.")
  }
  k_time <- seq(0, t_time, 1)
  tmp_sum <- sum(emgf(k_time, t_time, theta) * exp(
    theta * (t_time - k_time) * (-server_rate)))
  return((log(tmp_sum) - log(viol_prob - 1 + conflevel)) / theta)
}

#' optimizes the theta-dependent bound
#' @param blimit bandwidth limitation used for estimation of emgf.
statNC_optimize_ibl_theta <- function(viol_prob, t_time, conflevel, emgf,
                                      server_rate, blimit, accurate = FALSE) {

  step <- ifelse(accurate == FALSE, blimit / 1000, blimit / 10000)
  theta_vector <- seq(0 + step, blimit - step, step)
  backlog_vector <- sapply(theta_vector, statNC_inverse_backlog_theta,
                          viol_prob = viol_prob, t_time = t_time,
                          conflevel = conflevel, emgf = emgf,
                          server_rate = server_rate)
  # Sometimes the MGF gets too large -> Inf -> NaN -> Remove that
  # print(paste0("Stat: - which theta:",
  #             theta_vector[which.min(backlog_vector)],
  #             collapse = NULL))
  return(min(backlog_vector, na.rm = TRUE))
}

#' Compute an alternative interval for h_p
#' @param conflevel confidence level of hurst estimation
#' @param conflevel_beta compute another h_up for a higher confidence level.
statNC_interval <- function(
  sample_length, arrival_rate, hurst, std_dev, conflevel, iterations,
  viol_prob, t_time, server_rate, quantile_prob = 0.95,
  returnBLVector = FALSE) {

  # blimit is the bandwidth limitation of the estimator, in a single time step,
  # no more than blimit arrivals should occur.
  # 8 seems like a save choice, since the maximum in simulation is around ~5
  blimit <- 8
  iterations <- 10
  inverse_backlog <- rep(NA, iterations)

  for (i in 1:iterations) {
    f <- build_flow(
      arrival_rate = arrival_rate, hurst = hurst,
      sample_length = sample_length, std_dev = std_dev)

    emgf <- nonparametric_estimator(arrivals = f, conflevel = conflevel,
                                    blimit = blimit)
    inverse_backlog[i] <- statNC_optimize_ibl_theta(viol_prob = viol_prob,
                              t_time = t_time, conflevel = conflevel,
                              emgf = emgf, server_rate = server_rate,
                              blimit = blimit)

    .show_progress(i, iterations, "interval_h_up_quantile()")
  }

  if (returnBLVector) {
    return(inverse_backlog)
  } else {
    return(compute_statnc_interval(inverse_backlog, quantile_prob))
  }
}

compute_statnc_interval <- function(inverse_backlog, quantile_prob = 0.95) {
  ibl_mean <- mean(inverse_backlog)

  #cihelper = ci_help(inverse_backlog, quantile_prob)
  #lower = cihelper[[1]]
  #upper = cihelper[[2]]
  lower <- min(inverse_backlog)
  upper <- max(inverse_backlog)
  return(list("ibl_lower" = lower,
              "ibl_mean" = ibl_mean,
              "ibl_upper" = upper))
}

generate_values_csv <- function(
  sample_length, arrival_rate, hurst, time_n, server_rate, std_dev = 1.0,
  conflevel = 0.999, iterations = 10 ** 2) {

  d <- compute_distribution(
    arrival_rate = arrival_rate, hurst = hurst, sample_length = sample_length,
    time_n = time_n, server_rate = server_rate, std_dev = std_dev,
    iterations = iterations)

  # Take the "wrong" estimator now
  statnc_ibl <- statNC_interval(
    sample_length = sample_length, arrival_rate = arrival_rate,
    hurst = hurst, std_dev, conflevel = conflevel, iterations = iterations,
    viol_prob = 1 / iterations, t_time = time_n, server_rate = server_rate,
    returnBLVector = TRUE)

  df <- data.frame(bl_distribution = d, statnc_ibl)

  write.csv(df, file = "backlog_dist_statnc_fail.csv",
            col.names = TRUE, row.names = FALSE)
}

# Plots the empirical backlog distribution.

plot_distribution <- function(computed_dist, stat_mean, stat_lower, stat_upper,
                              stat_fail_mean, stat_fail_lower, stat_fail_upper,
                              conflevel, iterations, gran = 1000) {
  theme_set(theme_bw(base_size = 18))
  len <- length(computed_dist)
  maximum <- max(computed_dist)

  # Build the x-axis, start with 0 and end with the maximum
  bl <- seq(0, maximum, maximum / gran)
  # The cumulative backlog distribution curve
  # Init with 0
  pz <- rep(0, length(bl))
  labels  <-  data.frame(y = c(0.2, 0.4), x = c(stat_fail_mean, stat_mean),
                         label = c(round(stat_fail_mean, digits = 0),
                                   round(stat_mean, digits = 0)))

  # Build the cumulative distribution
  j  <-  1
  for (i in seq(0, maximum, maximum / gran)) {
    pz[j] <- length(computed_dist[computed_dist <= i]) / len
    j <- j + 1
  }
  # need violation probability, not confidence level
  nnb <- bl[min(which(pz >= 1 - (1 / iterations)))]

  frame <- data.frame(backlog = bl, perc = pz)
  # Prepare plot and plot backlog, trad and stat lines,
  # remove legend and set theme to bw
  q <- ggplot(frame, aes(x = backlog, y = perc)) +
    theme_bw(base_size = 18) +
    theme(legend.position = "none") +
    geom_line(size = 1, colour = "blue") +
    geom_vline(xintercept = c(nnb), colour = "blue") +
    # StatNC Fail Bound
    geom_vline(xintercept = c(stat_fail_mean), colour = "black") +
    geom_vline(xintercept = c(stat_fail_lower), colour = "black",
               linetype = "dotted") +
    geom_vline(xintercept = c(stat_fail_upper), colour = "black",
               linetype = "dotted") +
    # Correct StatNC Bound
    geom_vline(xintercept = c(stat_mean), colour = "aquamarine4") +
    geom_vline(xintercept = c(stat_lower), colour = "aquamarine4",
               linetype = "dotted") +
    geom_vline(xintercept = c(stat_upper), colour = "aquamarine4",
               linetype = "dotted") +
    geom_text(data = labels, aes(x = x, y = y, label = label)) +
    geom_label(aes(x = stat_fail_lower - 0.42 * maximum, y = 0.3,
                   label = "StatNC (IID)"), fill = "white", size = 5) +
    geom_label(aes(x = stat_fail_lower - 0.30 * maximum, y = 0.7,
                   label = "StatNC (LRD)"), fill = "white", size = 5) +
    geom_segment(aes(x = stat_fail_lower - 0.31 * maximum, y = 0.3,
                     xend = stat_fail_mean, yend = 0.3),
                 size = 0.4, arrow = NULL) +
    geom_segment(aes(x = stat_fail_lower - 0.03 * maximum, y = 0.7,
                     xend = stat_mean, yend = 0.7),
                 size = 0.4, arrow = NULL) +
    scale_x_log10() +
    annotation_logticks(sides = "b") +
    xlab("Backlog") +
    ylab("Cumulative Relative Frequencies")

  return(q)
}

# Computes the empirical backlog distribution and
# the corresponding traditional bound

plot_and_bound <- function(
  sample_length, arrival_rate, hurst, time_n, server_rate, std_dev = 1.0,
  splits = 20, conflevel = 0.999, iterations = 10 ** 2) {

  df <- read.csv(file = "backlog_dist_statnc_fail.csv", header = TRUE)

  snc_bound <- inverse_bound(
    time_n = time_n, std_dev = std_dev, hurst = hurst,
    arrival_rate = arrival_rate, server_rate = server_rate, p = 1 / iterations,
    splits = splits, conflevel = conflevel, estimated_h = FALSE)

  # Take the "wrong" estimator now
  # interval <- statNC_interval(
  # sample_length = sample_length, arrival_rate = arrival_rate,
  # hurst = hurst, std_dev, conflevel = conflevel, iterations = iterations,
  # viol_prob = 1 / iterations, t_time = time_n, server_rate = server_rate)

  fail_interval <- compute_statnc_interval(inverse_backlog = df$statnc_ibl)

  df_statnormal <- read.csv(file = "backlog_dist_h_confint.csv", header = TRUE)

  h.confint <- compute_h_up_quantile(h_vector = df_statnormal$hvector)

  stat_mean <- inverse_bound(
    time_n = time_n, std_dev = std_dev, hurst = h.confint$"Hurst_up_mean",
    arrival_rate = arrival_rate,
    server_rate = server_rate, p = 1 / iterations, splits = splits,
    conflevel = conflevel, estimated_h = TRUE)
  print(paste0("stat_mean = ", stat_mean))

  stat_lower <- inverse_bound(
    time_n = time_n, std_dev = std_dev, hurst = h.confint$"Hurst_lower_quant",
    arrival_rate = arrival_rate,
    server_rate = server_rate, p = 1 / iterations, splits = splits,
    conflevel = conflevel, estimated_h = TRUE)
  print(paste0("stat_lower = ", stat_lower))

  stat_upper <- inverse_bound(
    time_n = time_n, std_dev = std_dev, hurst = h.confint$"Hurst_upper_quant",
    arrival_rate = arrival_rate,
    server_rate = server_rate, p = 1 / iterations, splits = splits,
    conflevel = conflevel, estimated_h = TRUE)
  print(paste0("stat_upper = ", stat_upper))

  plot_distribution(
    computed_dist = df$bl_distribution, stat_mean = stat_mean,
    stat_lower = stat_lower,
    stat_upper = stat_upper,
    stat_fail_mean = fail_interval[["ibl_mean"]],
    stat_fail_lower = fail_interval[["ibl_lower"]],
    stat_fail_upper = fail_interval[["ibl_upper"]],
    conflevel = conflevel, iterations = iterations)
}

length_of_sample <- 2 ** 16
rate_arrival <- 10 ** (-2)
hurst_param <- 0.7
n_time <- 200
rate_server <- 1.5 * (10 ** (-2))
sigma_std <- 1.0
level_confidence <- 0.999
repetitions <- 500

# generate_values_csv(
#     sample_length = length_of_sample,
#     arrival_rate = rate_arrival, hurst = hurst_param, time_n = n_time,
#     server_rate = rate_server, std_dev = sigma_std,
#     conflevel = level_confidence, iterations = repetitions)

q <- plot_and_bound(
  sample_length = length_of_sample,
  arrival_rate = rate_arrival, hurst = hurst_param, time_n = n_time,
  server_rate = rate_server, std_dev = sigma_std, splits = 20,
  conflevel = level_confidence, iterations = repetitions)

pdf("backlog_distribution_StatNC_fail.pdf", width = 8, height = 5)

print(q)

dev.off()
