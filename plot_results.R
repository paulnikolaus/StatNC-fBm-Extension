##### plot_results.R #####
# File to visualize the results

# Plots the empirical backlog distribution.

plot_distribution <- function(computed_dist, stat, trad, gran = 1000) {
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
    geom_line(size = 1, colour = "#009E73") +
    geom_vline(xintercept = c(nnb)) +
    geom_vline(xintercept = c(trad), colour = "#56B4E9") +
    geom_vline(xintercept = c(stat), colour = "#E69F00") +
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
  arrival_rate, hurst, n, server_rate, std_dev = 1.0, splits = 20,
  conflevel = 0.995, iterations = 10 ** 2) {
  d <- compute_distribution(
    arrival_rate = arrival_rate, hurst = hurst, n = n,
    server_rate = server_rate, std_dev = std_dev, iterations = iterations)
  bound <- inverse_bound(
    n = n, std_dev = std_dev, hurst = hurst, arrival_rate = arrival_rate,
    server_rate = server_rate, p = 1 / iterations, splits = splits,
    conflevel = conflevel,
    traffic = NaN, estimate_traffic = FALSE)
  f <- build_flow(
    arrival_rate = arrival_rate, hurst = hurst, n = 2 ** 10, std_dev = std_dev)
  stat <- inverse_bound(
    n = n, std_dev = std_dev, hurst = hurst, arrival_rate = arrival_rate,
    server_rate = server_rate, p = 1 / iterations, splits = splits,
    conflevel = conflevel,
    traffic = f, estimate_traffic = TRUE)
  plot_distribution(d, stat, bound)
  
  #theme_set(theme_bw(base_size = 18))
  #qplot(x=1:length(d), y=d) + geom_line(aes(y=bound, colour="bound"))
  #return(list("SNC" = bound, "distribution" = d))
}

print(plot_and_bound(
  arrival_rate = 10 ** (-3), hurst = 0.7, n = 2*10 ** 2,
  server_rate = 5*10 ** (-3), std_dev = 1.0, splits = 20, conflevel = 0.999,
  iterations = 10 ** 3 - 1))
