# Load necessary library
library(dplyr)

# Read CSV files
quantMIS <- read.csv('Z:/segmentation/F8iia-quantification5.csv')
quantIM <- read.csv('Z:/segmentation/F8iic-quantificationV3.csv')

# MIS Ki67
# Calculate counts for each condition
CD3E_ki67 <- sum(quantMIS$CD3E > 250 & quantMIS$Ki67 > 160, na.rm = TRUE)
CD31_ki67 <- sum(quantMIS$CD31 > 578 & quantMIS$Ki67 > 160, na.rm = TRUE)
monocyte_ki67 <- sum(quantMIS$CD11c > 120 & quantMIS$CD11b < 300 & quantMIS$Ki67 > 160, na.rm = TRUE)
Ki67 <- sum(quantMIS$Ki67 > 160, na.rm = TRUE)

# Calculate percentages
CD3E_ki67_percentage <- (CD3E_ki67 / Ki67) * 100
CD31_ki67_percentage <- (CD31_ki67 / Ki67) * 100
monocyte_ki67_percentage <- (monocyte_ki67 / Ki67) * 100

# Print results
cat("CD3E Ki67 Percentage: ", CD3E_ki67_percentage, "%\n")
cat("CD31 Ki67 Percentage: ", CD31_ki67_percentage, "%\n")
cat("Monocyte Ki67 Percentage: ", monocyte_ki67_percentage, "%\n")

# IM Ki67
# Calculate counts for each condition
CD3E_ki67_IM <- sum(quantIM$CD3E > 240 & quantIM$Ki67 > 100, na.rm = TRUE)
tumor_ki67 <- sum(quantIM$MART1 > 800 & quantIM$Ki67 > 100, na.rm = TRUE)
Ki67_IM <- sum(quantIM$Ki67 > 100, na.rm = TRUE)
monocyte_ki67_IM <- sum(quantIM$CD11c > 100 & quantIM$CD11b < 220 & quantIM$Ki67 > 100, na.rm = TRUE)
macrophage_ki67 <- sum(quantIM$CD163 > 700 & quantIM$Ki67 > 100, na.rm = TRUE)
CD8_ki67 <- sum(quantIM$CD8a > 100 & quantIM$Ki67 > 100, na.rm = TRUE)

# Calculate percentages
tumor_ki67_percentage <- (tumor_ki67 / Ki67_IM) * 100
monocyte_ki67_percentage_IM <- (monocyte_ki67_IM / Ki67_IM) * 100
CD3E_ki67_percentage_IM <- (CD3E_ki67_IM / Ki67_IM) * 100
macrophage_ki67_percentage <- (macrophage_ki67 / Ki67_IM) * 100

# Print results
cat("Tumor Ki67 Percentage: ", tumor_ki67_percentage, "%\n")
cat("Monocyte Ki67 Percentage: ", monocyte_ki67_percentage_IM, "%\n")
cat("CD3E Ki67 Percentage: ", CD3E_ki67_percentage_IM, "%\n")
cat("Macrophage Ki67 Percentage: ", macrophage_ki67_percentage, "%\n")

library(ggplot2)
data_Ki67_stacked <- data.frame(
  category = c("MIS", "MIS", "MIS", "MIS","MIS","VGP","VGP", "VGP","VGP","VGP"),
  Subcategory = c("Endothelial cells", "Monocytes", "T cells", "Macrophage", "Tumor"),
  Percentage = c(2.727273, 28.18182, 33.63636, NA, NA, NA, 43.87384, 31.35867, 2.224019, 45.18803))

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
      axis.ticks = element_line()
  
    )
}

# Set color palette
MY_color=c("#d7191c","#fdae61", "#ffffbf", "#abd9e9","#2c7bb6")

plot <- ggplot(data_Ki67_stacked, aes(x = category, y = Percentage, fill = Subcategory)) +
  geom_bar(stat = "identity",position = position_dodge(width = 0.9))+
  geom_text(aes(label = sprintf("%.2f%%", Percentage)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3)+
  labs(title = "The cell type compositioin of Ki67+ cells in melanoma in-situ and invasive margin",
       x = "Region",
       y = "Percentage of Ki67+ per cell type")+ theme_jwmin() + scale_fill_manual(values=MY_color)

# Export an editable plot
save_pdf <- function(g, filename, width=10, height=10) {
  filepath <- paste0("~/Editable_plot/", filename)
  ggsave(plot=g, filename=filepath, width=width, height=height, bg="transparent",
         units="in", device=cairo_pdf)
}

save_pdf(plot, "The_cell_type_composition_of_Ki67+_cells.pdf")







