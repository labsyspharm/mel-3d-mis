library(patchwork)
library(ggplot2)

data_average_9um <- data.frame(
  Classification = c("9-18", "18-27", "27-35"),
  Percentage = c(96.80428,98.51973,96.41031))

data_average <- data.frame(
  Classification = c(9,18,27,35),
  Percentage = c(96.80428, 75.14495,39.81187,20.16689))

data_average_35um <- data.frame(
  Classification = c("18", "27", "35"),
  Percentage = c(75.14495,39.81187,20.16689))


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
MY_color=c("#fdae61","#ffffbf","#d7191c")
MY_color_1=c( "#e0f3f8","#91bfdb","#2c7bb6")
MY_color_2=c("#fdae61","#abd9e9","#2c7bb6","#d7191c")

# Plot
g1 <- ggplot(data_average_9um, aes(x = Classification, y = Percentage, fill = Classification)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = sprintf("%.2f%%", Percentage)), 
            position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  labs(x = "Z-range", y = "Percentage (%)") +
  scale_x_discrete(limits = c("9-18", "18-27", "27-35"))+ theme_jwmin() + scale_fill_manual(values=MY_color)


g2 <- ggplot(data_average_35um, aes(x = Classification, y = Percentage, fill = Classification)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = sprintf("%.2f%%", Percentage)), 
            position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  ylim(0, 100)+
  labs(x = NULL, y = "Percentage (%)") +
  scale_x_discrete(limits = c("9-27", "4-31","0-35"))+ theme_jwmin()+ scale_fill_manual(values=MY_color_1)


g3 <- ggplot(data_average, aes(x = Classification, y = Percentage, fill = Classification)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = sprintf("%.2f%%", Percentage)), 
            position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  geom_line(aes(x = as.numeric(Classification), y = -3.052 * as.numeric(Classification) + 125.880), 
            color = "#d7191c", lwd = 2) + 
  ylim(0, 100) +
  labs(x = "Z-range (μm)", y = "Percentage (%)", title = "Percentage of incomplete cells in synthetic 9 μm vs thicker sections") +
  scale_x_discrete(limits = c(9,18,27,35)) + theme_jwmin()  # Set limits if needed
  

# Combine plots with patchwork and add annotations
combined_plot <- g1 + g2 + 
  plot_layout(ncol = 2, width=c(15,15), axis_title="collect") + 
  plot_annotation(title = "Percentage of incomplete cells in synthetic 9 μm vs thick sections", 
                  theme = theme(plot.title = element_text(hjust = 0.5))) 
                  

# Export an editable plot
save_pdf <- function(g, filename, width=10, height=10) {
  filepath <- paste0("~/Editable_plot/", filename)
  ggsave(plot=g, filename=filepath, width=width, height=height, bg="transparent",
         units="in", device=cairo_pdf)
}

save_pdf(combined_plot, "Percentage_of_cells_missing_plot_with_more_bar.pdf")