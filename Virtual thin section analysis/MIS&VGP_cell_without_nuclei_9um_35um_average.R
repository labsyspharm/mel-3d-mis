library(patchwork)
library(ggplot2)

data_average <- data.frame(
  Classification = c("9","18","27", "35"),
  Percentage = c(12.36451,2.546152,1.105955, 0.7050216))

data_average <- data.frame(
  Classification = c("9", "35"),
  Percentage = c(12.36451, 0.7050216))

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
      axis.title.x = element_text(hjust = 0.5),
      axis.text.y = element_text(color = "black"),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(),
      legend.position="none"
    )
}

# Set color palette
MY_color=c( "#fdae61","#d7191c")


# Plot
g1 <- ggplot(data_average, aes(x = Classification, y = Percentage, fill = Classification)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = sprintf("%.2f%%", Percentage)), 
            position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  labs(x = "Z-range (μm)", y = "Percentage (%)", title = "Percentage of cells missing nuclei in synthetic 9 μm vs thicker sections") +
  scale_x_discrete(limits = c("9","35"))+ theme_jwmin() + scale_fill_manual(values=MY_color)
  

# Export an editable plot
save_pdf <- function(g, filename, width=10, height=10) {
  filepath <- paste0("~/Editable_plot/", filename)
  ggsave(plot=g, filename=filepath, width=width, height=height, bg="transparent",
         units="in", device=cairo_pdf)
}

save_pdf(g1, "Percentage_of_cells_missing_nuclei_plot-V1-YL.pdf")