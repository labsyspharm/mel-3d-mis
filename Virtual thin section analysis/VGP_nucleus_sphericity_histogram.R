# Load MIS and VGP datasets
VGP_35um_nuclei <- read.csv("Z:/F8ll_results/Invasive_margin_9_micron/melanoma-invasive_margin_nuclei_35um_nucleus_shape_(50vx)_stats.csv")
summary(VGP_35um_nuclei$height_length)
summary(VGP_35um_nuclei$height_width)


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
# VGP
library(ggplot2)
plot <- ggplot(VGP_35um_nuclei, aes(x = height_length)) + 
  geom_histogram(aes(y = after_stat(density)),  color = "black", fill = "white", binwidth = 0.05) +
  geom_vline(aes(xintercept = mean(height_length)), color = "royalblue", linetype = "dashed", linewidth = 1) +
  geom_density(alpha = 0.2, fill = "indianred") +
  labs(title = "The histogram of height-length ratio for intact nuclei in the melanoma invasive margin",
       x = "The height-length ratio (z/x) of nucleus bounding box in 35 micron", y = "Frequency") +
  scale_x_continuous(breaks = seq(0, 5, by=0.5))+ theme_jwmin() 

library(ggplot2)
plot <- ggplot(VGP_35um_nuclei, aes(x = height_width)) + 
  geom_histogram(aes(y = after_stat(density)),  color = "black", fill = "white", binwidth = 0.05) +
  geom_vline(aes(xintercept = mean(height_length)), color = "royalblue", linetype = "dashed", linewidth = 1) +
  geom_density(alpha = 0.2, fill = "indianred") +
  labs(title = "The histogram of height-width ratio for intact nuclei in melanoma invasive margin",
       x = "The height-width ratio (z/y) of nucleus bounding box in 35 micron", y = "Frequency") +
  scale_x_continuous(breaks = seq(0, 5, by=1))+ theme_jwmin() 

save_pdf <- function(g, filename, width=10, height=10) {
  filepath <- paste0("~/Editable_plot", filename)
  ggsave(plot=g, filename=filepath, width=width, height=height, bg="transparent",
         units="in", device=cairo_pdf)
}

save_pdf(plot, "The histogram of height-length ratio for intact nuclei in VGP.pdf")
