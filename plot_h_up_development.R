##### plot_h_up_development.R #####

library("reshape2") # melt
library("ggplot2")

plot_h_develop <- function(true_hurst = 0.7) {
  h_df <- read.csv(file = paste0(
    "h_up_developement_h_true=",
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

    geom_label(aes(
      x = 13000, y = max(h_df[, 2]) * 0.95,
      label = "Mean of H_up"
    ),
    fill = "white", size = 5
    ) +

    geom_label(aes(
      x = 5000, y = true_hurst * 0.99,
      label = "True Hurst Parameter"
    ),
    fill = "white", size = 5
    ) +

    theme_set(theme_bw(base_size = 19)) +
    theme(legend.position = "none") +
    xlab("Sample Size") +
    ylab(expression("Estimated H_up")) +
    theme(legend.title = element_blank())

  return(p)
}

ggsave("h_up_development.pdf", width = 8, height = 5, device = cairo_pdf)

print(plot_h_develop())

dev.off()
