##### plot_h_up_development.R #####

library("reshape2") # melt
library("ggplot2")

source("estimate_hurst.R")

# Show development of H_up for different sample sizes
#' @param true_hurst Value of the true hurst parameter
#' @return data frame of mean of h_up fir different sample sizes
h_development <- function(true_hurst = 0.7) {
  sample_sizes <- c(2**(10:17))
  h_ups <- rep(0.0, length(sample_sizes))
  for (i in 1:length(sample_sizes)) {
    h_vec <- est_h_up_vector(
      sample_length = sample_sizes[i], arrival_rate = 1.0,
      hurst = true_hurst, std_dev = 1.0, conflevel = 0.999,
      iterations = 500
    )
    h_ups[i] <- compute_h_up_quantile(h_vector = h_vec)$"Hurst_up_mean"
  }
  
  h_up_develop <- as.data.frame(cbind(sample_sizes, h_ups))
  
  write.csv(h_up_develop,
            file = paste0("h_up_development_h_true=", true_hurst, ".csv"),
            row.names = FALSE
  )
}

plot_h_develop <- function(true_hurst = 0.7) {
  h_df <- read.csv(file = paste0(
    "results/h_up_development_h_true=",
    true_hurst, ".csv"
  ))

  colnames(h_df) <- c(
    "SampleSize", "H_up"
  )

  long_df <- melt(h_df,
    id = "SampleSize",
    variable.name = "type",
    value.name = "H_up"
  )

  p <- ggplot(long_df, aes(
    x = SampleSize, y = H_up,
    group = type
  )) +
    scale_x_log10() +
    geom_line(aes(color = type, linetype = type), size = 0.8) +
    geom_point(aes(color = type, shape = type), size = 2.8) +
    scale_linetype_manual(values = "dotdash") +
    scale_color_manual(values = "blue") +
    scale_shape_manual(values = 20) +
    geom_hline(yintercept = true_hurst, linetype = "solid") +
    ylim(0.67, max(h_df[, 2])) +

    geom_label(aes(
      x = 14000, y = mean(h_df[, 2]),
      label = "Mean of H_up"
    ),
    fill = "white", size = 5
    ) +

    geom_label(aes(
      x = 5000, y = true_hurst * 0.99,
      label = "True Hurst Parameter H"
    ),
    fill = "white", size = 5
    ) +

    theme_set(theme_bw(base_size = 19)) +
    theme(legend.position = "none") +
    xlab("Sample Size") +
    ylab(expression("True H / Estimated H_up")) +
    theme(legend.title = element_blank())

  return(p)
}

print(h_development(true_hurst = 0.7))

# ggsave("results/h_up_development.pdf",
#   width = 8, height = 5,
#   device = cairo_pdf
# )
# 
# print(plot_h_develop(true_hurst = 0.7))
# 
# dev.off()
