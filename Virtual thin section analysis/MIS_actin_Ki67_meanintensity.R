library(ggplot2)
data_b_actin <- data.frame(
  Category = c("melanocytic cells", "T cells", "keratinocytes"),
  Percentage = c(624.2780613, 1190.212841, 423.5521486 ),
  sd =c(567.9825423, 704.0720782,497.4337263),
  sem = c(16.81481777, 15.86295194, 14.54260796))

data_Ki67 <- data.frame(
  Category = c("melanocytic cells", "T cells", "keratinocytes"),
  Percentage = c(28.46572033, 35.41805499,25.3937145),
  sd =c(36.73101161,86.81821683, 73.3753488),
  sem = c(1.087401849, 1.956040076, 2.145147936))

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
      legend.position = "None"
    )
}

# Set color palette
MY_color=c("#fc8d59", "#ffffbf", "#91bfdb")

# Create the bar graph
plot <- ggplot(data_b_actin, aes(x = Category, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin=Percentage-sem, ymax=Percentage+sem), width=.1,
                position=position_dodge(.9))+
  geom_text(aes(label = sprintf("%.2f", Percentage)), 
            position = position_dodge(width = 0.9), vjust = -1.5, size = 3) +
labs(title = "The mean intensity of b-actin for different cell types in the melanoma in-situ",
       x = "Cell type",
       y = "Mean intensity (GLU)") + theme_jwmin() + scale_fill_manual(values=MY_color)

# Export an editable plot
save_pdf <- function(g, filename, width=10, height=10) {
  filepath <- paste0("~/Editable_plot/", filename)
  ggsave(plot=g, filename=filepath, width=width, height=height, bg="transparent",
         units="in", device=cairo_pdf)
}

save_pdf(plot, "The_mean_intensity_of_b-actin_MIS.pdf")

plot1 <- ggplot(data_Ki67, aes(x = Category, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin=Percentage-sem, ymax=Percentage+sem), width=.1,
                position=position_dodge(.9))+
  geom_text(aes(label = sprintf("%.2f", Percentage)), 
            position = position_dodge(width = 0.9), vjust = -3.5, size = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), plot.title = element_text(hjust = 0.5)) +
  labs(title = "The mean intensity of Ki67 for different cell types in the melanoma in-situ",
       x = "Cell types",
       y = "Mean intensity (GLU)") + theme_jwmin() + scale_fill_manual(values=MY_color)

# Export an editable plot
save_pdf <- function(g, filename, width=10, height=10) {
  filepath <- paste0("~/Editable_plot/", filename)
  ggsave(plot=g, filename=filepath, width=width, height=height, bg="transparent",
         units="in", device=cairo_pdf)
}

save_pdf(plot1, "The_mean_intensity_of_Ki67_MIS.pdf")