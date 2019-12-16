##### plot_h_up_development.R #####

library("reshape2") # melt
library("ggplot2")

source("Bound.R")

# Show development of the backlog bound for different H_ups
#' @param true_hurst Value of the true hurst parameter
#' @return data frame of mean of h_up fir different sample sizes
backlog_development <- function(true_hurst = 0.7) {
  sample_h_ups <- read.csv(
    file = paste0("results/h_up_development_h_true=", true_hurst, ".csv"),
    header = T
  )
  h_ups <- sample_h_ups[, 2]
  backlog_bound_stat <- rep(0.0, length(h_ups))
  for (i in 1:length(h_ups)) {
    backlog_bound_stat[i] <- inverse_bound(
      time_n = 200, std_dev = 1.0, hurst = h_ups[i],
      arrival_rate = 10**(-2), server_rate = 1.5 * 10**(-2), p = 1 / 500,
      splits = 10, conflevel = 0.999,
      estimated_h = TRUE
    )
  }

  h_ups <- c(true_hurst, h_ups)
  backlog_bound_stat <- c(inverse_bound(
    time_n = 200, std_dev = 1.0, hurst = true_hurst,
    arrival_rate = 10**(-2), server_rate = 1.5 * 10**(-2), p = 1 / 500,
    splits = 10, conflevel = 0.999,
    estimated_h = FALSE
  ), backlog_bound_stat)

  backlog_bounds <- as.data.frame(cbind(h_ups, backlog_bound_stat))

  write.csv(backlog_bounds,
    file = paste0("results/backlog_development_h_true=", true_hurst, ".csv"),
    row.names = FALSE
  )
}

plot_backlog_develop <- function(true_hurst = 0.7) {
  backlog_df <- read.csv(file = paste0(
    "results/backlog_development_h_true=",
    true_hurst, ".csv"
  ))

  true_backlog <- backlog_df[1, 2]
  backlog_df <- backlog_df[-1, ]

  colnames(backlog_df) <- c(
    "H_up", "BacklogBound"
  )

  long_df <- melt(backlog_df,
    id = "H_up",
    variable.name = "type",
    value.name = "BacklogBound"
  )

  p <- ggplot(long_df, aes(
    x = H_up, y = BacklogBound,
    group = type
  )) +
    geom_line(aes(color = type, linetype = type), size = 0.8) +
    geom_point(aes(color = type, shape = type), size = 2.8) +
    scale_x_reverse() +
    scale_linetype_manual(values = "dotdash") +
    scale_color_manual(values = "blue") +
    scale_shape_manual(values = 20) +
    # ylim(0.67, 0.8) +

    geom_hline(yintercept = true_backlog, linetype = "solid") +
    geom_label(aes(
      x = 0.76, y = mean(backlog_df[, 2]),
      label = "StatNC Backlog Bound"
    ),
    fill = "white", size = 5
    ) +

    geom_label(aes(
      x = 0.75, y = true_backlog * 0.95,
      label = "SNC Backlog Bound"
    ),
    fill = "white", size = 5
    ) +

    theme_bw(base_size = 19) +
    theme(legend.position = "none") +
    xlab("Estimated Hurst Parameter") +
    ylab("Backlog Bound") +
    theme(legend.title = element_blank())

  return(p)
}

# print(backlog_development(true_hurst = 0.7))

ggsave("results/backlog_development.pdf",
  width = 8, height = 5,
  device = cairo_pdf
)

print(plot_backlog_develop(true_hurst = 0.7))

dev.off()
