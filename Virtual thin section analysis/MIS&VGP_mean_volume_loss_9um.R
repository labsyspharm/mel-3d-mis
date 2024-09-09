library(ggplot2)

# Set theme
theme_jwmin <- function() {
  theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black"),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA) ,
      plot.title = element_text(hjust = 0.5),
      text = element_text(family = "Arial"),
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(),
      legend.position="none"
    )
}

# Set color palette
MY_color=c("#fc8d59", "#ffffbf", "#91bfdb")

data_mean_volume_loss_avg <- data.frame(
  Classification = c("All cells", "CD3+ T cells", "MART1+ tumor cells"),
  Percentage = c(50.00813223, 47.06484973, 47.8587892),
  IQR = c(58.64925614,58.45818449, 56.26461137))


plot <- ggplot(data_mean_volume_loss_avg, aes(x = Classification, y = Percentage, fill = Classification)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin=Percentage-IQR/2, ymax=Percentage+IQR/2), width=.1,
                position=position_dodge(.9))+
  geom_text(aes(label = sprintf("%.2f%%", Percentage)), 
            position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  labs(title = "Average cell volume missing from synthetic 9 μm sections",
       x = "9 μm section",
       y = "Percentage (%)") + theme_jwmin() + scale_fill_manual(values=MY_color)

# Export an editable plot
save_pdf <- function(g, filename, width=10, height=10) {
  filepath <- paste0("~/Editable_plot/", filename)
  ggsave(plot=g, filename=filepath, width=width, height=height, bg="transparent",
         units="in", device=cairo_pdf)
}

save_pdf(plot, "Average_cell_volume_missing_plot.pdf")

