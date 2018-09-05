##### plot_utilization.R #####

library("reshape2") # melt
library("ggplot2")

source("Bound.R") # inverse_bound(), loads estimate_hurst.R and simulation.R

# Plots the the bound against the utilization.

csv_backlog_vs_util <- function(
  sample_length, arrival_rate, hurst, time_n, conflevel, prob, iterations,
  std_dev = 1.0, splits = 20) {
  # NOTE: all 1 / iterations have been replaced with prob
  utilizations <- (14:19) / 20

  if ((prob) < (1 - conflevel)) {
    stop(paste0(
      "p = ", prob, " < (1 - conflevel) = ",
      1 - conflevel, ". \n
    The bound runs in an infinite loop as the stat_backlog_bound() bound can
    never be below (1-alpha)"
    ))
  }

  simulated_backlog <- rep(NA, length(utilizations))
  snc_bound <- rep(NA, length(utilizations))
  stat_mean <- rep(NA, length(utilizations))
  stat_low <- rep(NA, length(utilizations))
  stat_up <- rep(NA, length(utilizations))
  i <- 1

  for (util in utilizations) {
    print(paste0("utilization: ", util))
    simulated_backlog[i] <- quantile((compute_distribution(
      arrival_rate = arrival_rate, hurst = hurst, sample_length = sample_length,
      time_n = time_n, server_rate = 1 / util, std_dev = std_dev,
      iterations = iterations
    )), probs = 1 - (prob))
    print(paste0("simulated_backlog: ", simulated_backlog[i]))

    snc_bound[i] <- inverse_bound(
      time_n = time_n, std_dev = std_dev, hurst = hurst,
      arrival_rate = arrival_rate, server_rate = 1 / util, p = prob,
      splits = splits, conflevel = conflevel, estimated_h = FALSE
    )
    print(paste0("snc_bound: ", snc_bound[i]))

    h_up_vec <- est_h_up_vector(
      sample_length = sample_length, arrival_rate = arrival_rate,
      hurst = hurst, std_dev = std_dev, conflevel = conflevel,
      iterations = iterations
    )

    h_up_quantile <- compute_h_up_quantile(
      h_vector = h_up_vec,
      quantile_prob = 0.95
    )
    print(h_up_quantile)

    stat_mean[i] <- inverse_bound(
      time_n = time_n, std_dev = std_dev,
      hurst = h_up_quantile$"Hurst_up_mean",
      arrival_rate = arrival_rate, server_rate = 1 / util, p = prob,
      splits = splits, conflevel = conflevel, estimated_h = TRUE
    )
    print(paste0("stat_mean: ", stat_mean[i]))

    stat_low[i] <- inverse_bound(
      time_n = time_n, std_dev = std_dev,
      hurst = h_up_quantile$"Hurst_lower_quant",
      arrival_rate = arrival_rate, server_rate = 1 / util, p = prob,
      splits = splits, conflevel = conflevel, estimated_h = TRUE
    )

    stat_up[i] <- inverse_bound(
      time_n = time_n, std_dev = std_dev,
      hurst = h_up_quantile$"Hurst_upper_quant",
      arrival_rate = arrival_rate, server_rate = 1 / util, p = prob,
      splits = splits, conflevel = conflevel, estimated_h = TRUE
    )

    i <- i + 1
  }


  backlog_bounds_df <- as.data.frame(
    cbind(
      utilizations, stat_up, stat_mean, stat_low, snc_bound,
      simulated_backlog
    )
  )

  write.csv(backlog_bounds_df,
    file = "backlog_bounds.csv",
    row.names = FALSE
  )

  return(backlog_bounds_df)
}

plot_backlog_vs_util <- function() {
  backlog_bounds_df <- read.csv(file = "backlog_bounds.csv")

  colnames(backlog_bounds_df) <- c(
    "utilizations", "StatNC up", "Mean of StatNC bounds", "StatNC low",
    "SNC Bound", "Simulation"
  )

  long_df <- melt(backlog_bounds_df,
    id = "utilizations",
    variable.name = "type",
    value.name = "Backlog_bound"
  )

  p <- ggplot(long_df, aes(
    x = utilizations, y = Backlog_bound,
    group = type
  )) +
    geom_line(aes(color = type, linetype = type), size = 0.8) +
    geom_point(aes(color = type, shape = type), size = 2.8) +
    scale_linetype_manual(
      values = c("dashed", "solid", "dashed", "F1", "dotdash")
    ) +
    scale_color_manual(
      values = c("aquamarine4", "black", "aquamarine4", "red", "blue")
    ) +
    scale_shape_manual(values = c(20, 19, 20, 18, 17)) +
    ylim(0.5, max(backlog_bounds_df)) +

    geom_label(aes(
      x = 0.83, y = max(backlog_bounds_df[3]) * 0.85,
      label = "Mean of StatNC bounds"
    ),
    fill = "white", size = 5
    ) +
    geom_label(aes(
      x = 0.85, max(backlog_bounds_df[5]) * 0.6,
      label = "SNC Bound"
    ),
    fill = "white", size = 5
    ) +
    geom_label(aes(
      x = 0.93, max(backlog_bounds_df[6]) * 0.6,
      label = "Simulation"
    ), fill = "white", size = 5) +

    theme_set(theme_bw(base_size = 19)) +
    # theme(legend.position = c(0.25, 0.8),
    #       legend.background = element_rect(color = "black"),
    #       axis.text = element_text(size = 20)) +
    theme(legend.position = "none") +
    xlab("Utilization") +
    ylab("Backlog") +
    theme(legend.title = element_blank())

  return(p)
}

# csv_backlog_vs_util(
#   sample_length = 2 ** 16,
#   arrival_rate = 10 ** (-2), hurst = 0.7, time_n = 200, conflevel = 0.999,
#   prob = 1 / 500, iterations = 5000, std_dev = 1.0, splits = 20)

# pdf("backlog_vs_util.pdf", width = 8, height = 5)
ggsave("backlog_vs_util.pdf", width = 8, height = 5, device = cairo_pdf)

plot_backlog_vs_util()

dev.off()
