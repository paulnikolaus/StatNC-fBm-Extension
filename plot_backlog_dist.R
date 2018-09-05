##### plot_backlog_dist.R #####

library("ggplot2")

source("simulation.R") #' loads compute_distribution()
source("estimate_hurst.R") # loads the necessary tools for estimation
source("Bound.R") #' loads inverse_bound()

generate_values_csv <- function(
  sample_length, arrival_rate, hurst, time_n, server_rate, std_dev = 1.0,
  conflevel = 0.999, iterations = 10**2) {
  d <- compute_distribution(
    arrival_rate = arrival_rate, hurst = hurst, sample_length = sample_length,
    time_n = time_n, server_rate = server_rate, std_dev = std_dev,
    iterations = iterations
  )

  hvector <- est_h_up_vector(
    sample_length = sample_length, arrival_rate = arrival_rate, hurst = hurst,
    std_dev = std_dev, conflevel = conflevel, iterations = iterations
  )

  df <- data.frame(bl_distribution = d, hvector = hvector)

  write.csv(df,
    file = "results/backlog_dist_h_confint.csv",
    col.names = TRUE, row.names = FALSE
  )
}

# Plots the empirical backlog distribution.

plot_distribution <- function(computed_dist, stat_mean, stat_lower, stat_upper,
                              trad, conflevel, iterations, gran = 1000) {
  theme_set(theme_bw(base_size = 18))
  len <- length(computed_dist)
  maximum <- max(computed_dist)

  # Build the x axis, start with 0 and end with the maximum
  bl <- seq(0, maximum, maximum / gran)
  # The cumulative backlog distribution curve
  # Init with 0
  pz <- rep(0, length(bl))
  labels <- data.frame(
    y = c(0.2, 0.4), x = c(trad, stat_mean),
    label = c(
      round(trad, digits = 0),
      round(stat_mean, digits = 0)
    )
  )

  # Build the cumulative distribution
  j <- 1
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
    geom_vline(xintercept = c(trad), colour = "red") +
    geom_vline(xintercept = c(stat_mean), colour = "black") +
    geom_vline(
      xintercept = c(stat_lower), colour = "aquamarine4",
      linetype = "dotted"
    ) +
    geom_vline(
      xintercept = c(stat_upper), colour = "aquamarine4",
      linetype = "dotted"
    ) +
    geom_text(data = labels, aes(x = x, y = y, label = label)) +
    geom_label(aes(x = nnb - 0.4 * maximum, y = 0.3, label = "SNC"),
      fill = "white", size = 5
    ) +
    geom_label(aes(x = nnb - 0.4 * maximum, y = 0.7, label = "StatNC"),
      fill = "white", size = 5
    ) +
    # annotate("text", x = c(nnb - 0.4 * maximum, nnb - 0.4 * maximum),
    #         y = c(0.30, 0.70), label = c("SNC", "StatNC"), size = 5.5) +
    geom_segment(aes(
      x = nnb - maximum / 7, y = 0.3, xend = trad,
      yend = 0.3
    ), size = 0.4, arrow = NULL) +
    geom_segment(aes(
      x = nnb - maximum / 12, y = 0.7, xend = stat_mean,
      yend = 0.7
    ), size = 0.4, arrow = NULL) +
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
                           splits = 20, conflevel = 0.999, iterations = 10**2) {
  if ((1 / iterations) < (1 - conflevel)) {
    stop(paste0(
      "p = ", (1 / iterations), " < (1 - conflevel) = ",
      1 - conflevel, ". \n
      The bound runs in an infinite loop as the stat_backlog_bound() bound can
      never be below (1-alpha)"
    ))
  }

  df <- read.csv(file = "results/backlog_dist_h_confint.csv", header = T)

  h.confint <- compute_h_up_quantile(h_vector = df$hvector)

  snc_bound <- inverse_bound(
    time_n = time_n, std_dev = std_dev, hurst = hurst,
    arrival_rate = arrival_rate, server_rate = server_rate, p = 1 / iterations,
    splits = splits, conflevel = conflevel, estimated_h = FALSE
  )

  print(paste0(
    "Hurst_lower_quant = ", h.confint[1],
    ", Hurst_up_mean = ", h.confint[2],
    ", Hurst_upper_quant = ", h.confint[3]
  ))

  stat_mean <- inverse_bound(
    time_n = time_n, std_dev = std_dev, hurst = h.confint$"Hurst_up_mean",
    arrival_rate = arrival_rate,
    server_rate = server_rate, p = 1 / iterations, splits = splits,
    conflevel = conflevel, estimated_h = TRUE
  )
  print(paste0("stat_mean = ", stat_mean))

  stat_lower <- inverse_bound(
    time_n = time_n, std_dev = std_dev, hurst = h.confint$"Hurst_lower_quant",
    arrival_rate = arrival_rate,
    server_rate = server_rate, p = 1 / iterations, splits = splits,
    conflevel = conflevel, estimated_h = TRUE
  )
  print(paste0("stat_lower = ", stat_lower))

  stat_upper <- inverse_bound(
    time_n = time_n, std_dev = std_dev, hurst = h.confint$"Hurst_upper_quant",
    arrival_rate = arrival_rate,
    server_rate = server_rate, p = 1 / iterations, splits = splits,
    conflevel = conflevel, estimated_h = TRUE
  )
  print(paste0("stat_upper = ", stat_upper))

  plot_distribution(
    computed_dist = df$bl_distribution, stat_mean = stat_mean,
    stat_lower = stat_lower, stat_upper = stat_upper, trad = snc_bound,
    conflevel = conflevel, iterations = iterations
  )

  #' theme_set(theme_bw(base_size = 18))
  #' qplot(x = 1:length(d), y = d) +
  #' geom_line(aes(y = bound, color = "bound"))
  #' return(list("SNC" = bound, "distribution" = d))
}

length_of_sample <- 2**16
rate_arrival <- 10**(-2)
hurst_param <- 0.7
n_time <- 200
rate_server <- 1.5 * (10**(-2))
sigma_std <- 1.0
repetitions <- 500
level_confidence <- 0.999

# generate_values_csv(sample_length = length_of_sample,
#   arrival_rate = rate_arrival, hurst = hurst_param, time_n = n_time,
#   server_rate = 1.5 * (10 ** (-2)), std_dev = sigma_std,
#   conflevel = level_confidence, iterations = repetitions)

q <- plot_and_bound(
  sample_length = length_of_sample,
  arrival_rate = rate_arrival, hurst = hurst_param, time_n = n_time,
  server_rate = 1.5 * (10**(-2)), std_dev = sigma_std, splits = 20,
  conflevel = level_confidence, iterations = repetitions
)

# pdf("results/backlog_distribution.pdf", width = 8, height = 5)
ggsave("results/backlog_distribution.pdf", width = 8, height = 5,
  device = cairo_pdf)

print(q)

dev.off()
