##### plot_backlog_dist.R #####
# File to visualize the results

library("ggplot2")

source("simulation.R") # compute_distribution(), loads bound.R
source("estimate_hurst.R") # loads the necessary tools for estimation
source("Bound.R") # inverse_bound()

# Plots the empirical backlog distribution.

plot_distribution <- function(computed_dist, stat, stat_lower, stat_upper,
                              trad, gran = 1000) {
  theme_set(theme_bw(base_size = 18))
  len <- length(computed_dist)
  maximum <- max(computed_dist)

  # Build the x axis, start with 0 and end with the maximum
  bl <- seq(0, maximum, maximum / gran)
  # The cumulative backlog distribution curve
  # Init with 0
  pz <- rep(0, length(bl))
  labels  <-  data.frame(y = c(0.2, 0.4), x = c(trad, stat),
                         label = c(round(trad, digits = 3),
                                   round(stat, digits = 3)))

  # Build the cumulative distribution
  j  <-  1
  for (i in seq(0, maximum, maximum / gran)) {
    pz[j] <- length(computed_dist[computed_dist <= i]) / len
    j <- j + 1
  }
  # 99,99% percentile
  nnb <- bl[min(which(pz >= 0.9999))]

  frame <- data.frame(backlog = bl, perc = pz)
  # Prepare plot and plot backlog, trad and stat lines,
  # remove legend and set theme to bw

  q <- ggplot(frame, aes(x = backlog, y = perc)) +
    theme_bw(base_size = 18) +
    theme(legend.position = "none") +
    geom_line(size = 1, colour = "blue") +
    geom_vline(xintercept = c(nnb), colour = "blue") +
    geom_vline(xintercept = c(trad), colour = "red") +
    geom_vline(xintercept = c(stat), colour = "black") +
    geom_vline(xintercept = c(stat_lower), colour = "darkgreen",
               linetype = "dotted") +
    geom_vline(xintercept = c(stat_upper), colour = "#lightgreen",
               linetype = "dotted") +
    geom_text(data = labels, aes(x = x, y = y, label = label)) +
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
  splits = 20, conflevel = 0.995, iterations = 10 ** 2) {
  d <- compute_distribution(
    arrival_rate = arrival_rate, hurst = hurst, sample_length = sample_length,
    time_n = time_n, server_rate = server_rate, std_dev = std_dev,
    iterations = iterations)
  bound <- inverse_bound(
    time_n = time_n, std_dev = std_dev, hurst = hurst,
    arrival_rate = arrival_rate, server_rate = server_rate, p = 1 / iterations,
    splits = splits, conflevel = conflevel, estimated_h = FALSE)


  # h.confint <- confint_of_h_up(
  #   sample_length = sample_length, arrival_rate = arrival_rate, hurst = hurst,
  #   std_dev = std_dev, conflevel = conflevel, iterations = iterations,
  #   confint.conflevel = 0.95)

  # c(h_estimated, h_up, h_up^beta) from interval_h_up_alter()
  h.confint <- interval_h_up_alter(
    sample_length = sample_length, arrival_rate = arrival_rate, hurst = hurst,
    std_dev = std_dev, conflevel = conflevel, iterations = iterations,
    conflevel_beta = 0.99999)

  # h_up <- flow_to_h_up(f, arrival_rate = arrival_rate, std_dev = std_dev,
  #                      conflevel = conflevel)
  # print(paste0("Hurst_up_mean = ", h.confint[1],
  #              ", Hurst_up_lower = ", h.confint[2],
  #              ", Hurst_up_upper = ", h.confint[3]))
  print(paste0("Hurst_mean = ", h.confint[1],
               ", Hurst_up_mean = ", h.confint[2],
               ", Hurst_up_beta = ", h.confint[3]))

  stat_mean <- inverse_bound(
    time_n = time_n, std_dev = std_dev, hurst = h.confint$"hurst_up",
    arrival_rate = arrival_rate,
    server_rate = server_rate, p = 1 / iterations, splits = splits,
    conflevel = conflevel, estimated_h = TRUE)
  print(paste0("stat_mean = ", stat_mean))

  stat_lower <- inverse_bound(
    time_n = time_n, std_dev = std_dev, hurst = h.confint$"hurst_est",
    arrival_rate = arrival_rate,
    server_rate = server_rate, p = 1 / iterations, splits = splits,
    conflevel = conflevel, estimated_h = TRUE)
  print(paste0("stat_lower = ", stat_lower))

  stat_upper <- inverse_bound(
    time_n = time_n, std_dev = std_dev, hurst = h.confint$"hurst_up_beta",
    arrival_rate = arrival_rate,
    server_rate = server_rate, p = 1 / iterations, splits = splits,
    conflevel = conflevel, estimated_h = TRUE)
  print(paste0("stat_upper = ", stat_upper))

  plot_distribution(
    computed_dist = d, stat = stat_mean, stat_lower = stat_lower,
    stat_upper = stat_upper, trad = bound)

  # theme_set(theme_bw(base_size = 18))
  # qplot(x = 1:length(d), y = d) +
  # geom_line(aes(y = bound, color = "bound"))
  # return(list("SNC" = bound, "distribution" = d))
}

q <- plot_and_bound(
  sample_length = 2 ** 15,
  arrival_rate = (10 ** (-2)), hurst = 0.7, time_n = 200,
  server_rate = 1.5 * (10 ** (-2)), std_dev = 1.0, splits = 20,
  conflevel = 0.999, iterations = 200)
pdf("backlog_distribution.pdf", width = 8, height = 5)

print(q)

# results:
# blue line (SNC-bound): 164.0
# yellow line (StatNC-bound): 204.5

dev.off()
