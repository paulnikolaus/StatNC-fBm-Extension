##### plot_utilization.R #####
# File to visualize the results

library("reshape2")  # melt
library("ggplot2")

source("simulation.R") # compute_distribution(), loads bound.R
source("estimate_hurst.R") # loads the necessary tools for estimation
source("Bound.R") # inverse_bound()

# Plots the the bound against the utilization.

backlog_vs_util <- function(
  sample_length, arrival_rate, hurst, time_n, std_dev = 1.0,
  splits = 20, conflevel = 0.995, iterations = 10 ** 2) {
  utilizations <- (10:19) / 20

  simulated_backlog <- rep(NA, length(utilizations))
  bound <- rep(NA, length(utilizations))
  h_up_est <- rep(NA, length(utilizations))
  stat_mean <- rep(NA, length(utilizations))
  i <- 1

  for (util in utilizations) {
    print(paste0("utilization: ", util))
    simulated_backlog[i] <- quantile((compute_distribution(
      arrival_rate = arrival_rate, hurst = hurst, sample_length = sample_length,
      time_n = time_n, server_rate = 1 / util, std_dev = std_dev,
      iterations = iterations)), probs = conflevel)
    print(paste0("simulated_backlog: ", simulated_backlog[i]))

    bound[i] <- inverse_bound(
      time_n = time_n, std_dev = std_dev, hurst = hurst,
      arrival_rate = arrival_rate, server_rate = 1 / util, p = 1 / iterations,
      splits = splits, conflevel = conflevel, estimated_h = FALSE)
    print(paste0("bound: ", bound[i]))

    h_up_est[i] <- mean_of_h_up(
      sample_length = sample_length, arrival_rate = arrival_rate, hurst = hurst,
      std_dev = std_dev, conflevel = conflevel, iterations = iterations)

    stat_mean[i] <- inverse_bound(
      time_n = time_n, std_dev = std_dev, hurst = h_up_est[i],
      arrival_rate = arrival_rate,
      server_rate = 1 / util, p = 1 / iterations, splits = splits,
      conflevel = conflevel, estimated_h = TRUE)
    print(paste0("stat_mean: ", stat_mean[i]))

    i <- i + 1
  }


  backlog_bounds_df = as.data.frame(
    cbind(utilizations, stat_mean, bound, simulated_backlog))

  colnames(backlog_bounds_df) <- c(
    "utilizations", "Mean of StatNC bounds", "SNC Bound",
    "Quantile of Simulated Backlog")

  # print(backlog_bounds_df)

  long_df <- melt(backlog_bounds_df, id = "utilizations",
                  variable.name = "type",
                  value.name = "Backlog_bound")
  # print(long_df)

  p <- ggplot(long_df, aes(x = utilizations, y = Backlog_bound,
                                  group = type)) +
    geom_line(aes(color = type, linetype = type), size = 0.8) +
    geom_point(aes(color = type, shape = type), size = 2.8) +
    scale_linetype_manual(values = c("dashed", "F1", "dotdash")) +
    scale_color_manual(values = c("black", "red", "blue")) +
    scale_shape_manual(values = c(17, 19, 18)) +
    theme_set(theme_bw(base_size = 19)) +
    theme(legend.position = c(0.25, 0.8),
          legend.background = element_rect(color = "black"),
          axis.text = element_text(size = 18)) +
    xlab("Utilization") +
    ylab("Backlog Bound") +
    theme(legend.title = element_blank())

  return(p)
}

pdf("backlog_vs_util.pdf", width = 8, height = 5)

backlog_vs_util(
  sample_length = 2 ** 14,
  arrival_rate = 10 ** (-3), hurst = 0.7, time_n = 2 * (10 ** 2),
  std_dev = 1.0, splits = 20, conflevel = 0.999,
  iterations = 2 * (10 ** 2))

dev.off()
