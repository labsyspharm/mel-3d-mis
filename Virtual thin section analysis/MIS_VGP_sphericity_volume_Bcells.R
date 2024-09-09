# Load necessary libraries
library(ggplot2)

# Sphericity data
sphericity <- data.frame(
  Condition = c("MIS", "VGP"),
  sphericity = c(0.347666667, 0.583639708),
  std_sphericity= c(0.065565913, 0.10169758))

# Volume data
volume <- data.frame(
  Condition = c("MIS", "VGP"),
  volume = c(316.9212667, 160.9151947),
  std_volume= c(192.9894597,143.1328829))

# Set theme
theme_jwmin <- function() {
  theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black"),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA) ,
      text = element_text(family = "Arial"),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(),
      legend.position="none"
    )
}


plot1 <- ggplot(sphericity, aes(x = Condition, y = sphericity)) +
  geom_bar(stat = "identity", fill = "#4393c3", width = 0.5) +
  geom_errorbar(aes(ymin = sphericity - std_sphericity, ymax = sphericity + std_sphericity), width = 0.2) +
  geom_text(aes(label = round(sphericity, 2)), 
            position = position_dodge(width=0.9), vjust = -0.5, size = 5)+
  ggtitle("The sphericity of B cells in the MIS and VGP") +
  labs(x = "Region", y ="Mean Sphericity") + theme_jwmin() 


# Plot Volume
plot2 <- ggplot(volume, aes(x = Condition, y = volume)) +
  geom_bar(stat = "identity", fill = "#d6604d", width = 0.5) +
  geom_errorbar(aes(ymin = volume - std_volume, ymax = volume + std_volume), width = 0.2) +
  geom_text(aes(label = round(volume,2)), 
            position = position_dodge(width=0.9),  vjust = -0.5, size = 5)+
  ggtitle("The volume of B cells in the MIS and VGP") +
  labs(x = "Region", y = expression("Mean volume" ~(μm^3))) + theme_jwmin() 

# Export an editable plot
save_pdf <- function(g, filename, width=10, height=10) {
  filepath <- paste0("~/Editable_plot/", filename)
  ggsave(plot=g, filename=filepath, width=width, height=height, bg="transparent",
         units="in", device=cairo_pdf)
}

save_pdf(plot1, "The_sphericity_B_cells_plot.pdf")
save_pdf(plot2, "The_volume_B_cells_plot.pdf")