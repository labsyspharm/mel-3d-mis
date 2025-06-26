# Load necessary libraries
library(ggplot2)
library(dplyr)
MHC1_MX1 <- read.csv("MHC1vsMX1.csv")
upper_bound <- quantile(MHC1_MX1$MHC1_intensity, 0.95, na.rm = TRUE)
lower_bound <- quantile(MHC1_MX1$MHC1_intensity, 0.05, na.rm = TRUE)
# Filter out the 1% outliers
MHC1_MX1_filtered <- MHC1_MX1 %>%
  filter(MHC1_intensity>= lower_bound & MHC1_intensity<= upper_bound)


# Summarize stats
MHC1_stats <- MHC1_MX1_filtered %>%
  group_by(Condition) %>%
  summarise(mean = mean(MHC1_intensity, na.rm = TRUE),
            sd = sd(MHC1_intensity, na.rm = TRUE),
            n = n())

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

plot <- ggplot(MHC1_stats , aes(x = Condition, y = mean/100)) +
  geom_bar(stat = "identity", fill = "#4393c3", width = 0.5) +
  geom_jitter(data = MHC1_MX1_filtered , aes(x = Condition, y = MHC1_intensity/100), 
              width = 0.1, size = 1, color = "gray") +
  geom_errorbar(aes(ymin = (mean - sd)/100, ymax = (mean + sd)/100, width = 0.1)) +
  # geom_text(aes(label = round(mean, 2)), 
            # position = position_dodge(width=0.9), vjust = -0.5, size = 5)+
  ggtitle("MHC-1 expression in relation to MX1 puncta") +
  labs(x = "Proximity to MX1 puncta", y = expression("MHC-1 Mean Intensity" ~times~ 10^2 ~ "(GLU)")) + theme_jwmin() 


# Export an editable plot
save_pdf <- function(g, filename, width=10, height=10) {
  filepath <- paste0("~/Editable_plot/", filename)
  ggsave(plot=g, filename=filepath, width=width, height=height, bg="transparent",
         units="in", device=cairo_pdf)
}

save_pdf(plot, "MHC-1 expression in relation to MX1 puncta_plot.pdf")
