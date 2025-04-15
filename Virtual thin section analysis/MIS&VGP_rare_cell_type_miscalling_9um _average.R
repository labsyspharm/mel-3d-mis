library(ggplot2)
# 
# # MIS and VGP
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
MY_color=c( "#fee090","#abd9e9","#e0f3f8","#fdae61","#d73027","#74add1","#f46d43","#4575b4")


# Plot data
# MIS and VGP
data_miscall <- data.frame(
  Classification = c("LAG3+ cells", "MX1+ cells", "LAG3+ & MX1+ cells","CD103+ cells","CD8+ cells","CD4+ cells","MART1+ cells", "SOX10+ cells"),
  Percentage = c(36.5775, 36.57916,23.83259, 13.52474,12.69739, 8.730851, 5.328166, 1.343655))

plot <- ggplot(data_miscall, aes(x = Classification, y = Percentage, fill = Classification)) +
  geom_bar(stat = "identity")+
  geom_text(aes(label = sprintf("%.2f%%", Percentage)), position = position_dodge(width = 0.9),vjust = -0.5)+
  labs(title = "Percentage of misclassified cell types in thin sections",
       x = "Cell type",
       y = "Percentage (%)")+
scale_x_discrete(limits = c("LAG3+ cells", "MX1+ cells", "LAG3+ & MX1+ cells","CD103+ cells","CD8+ cells","CD4+ cells","MART1+ cells", "SOX10+ cells")) + theme_jwmin() + scale_fill_manual(values=MY_color)

# Export an editable plot
save_pdf <- function(g, filename, width=10, height=10) {
  filepath <- paste0("~/Editable_plot", filename)
  ggsave(plot=g, filename=filepath, width=width, height=height, bg="transparent",
         units="in", device=cairo_pdf)
}

save_pdf(plot, "Percentage_of_misclassified_cell_types_plot_V2.pdf")