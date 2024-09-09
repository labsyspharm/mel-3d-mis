library(dplyr)

# Read CSV files
distransform_MIS <- read.csv('Z:/F8ll_results/MIS_distance_transform.csv')
distransform_VGP <- read.csv('Z:/F8ll_results/VGP_distance_transform.csv')
# MIS 
# Calculate counts for each condition
MIS_incomplete_cell <- sum(distransform_MIS$MeanIntensity <= 1.5 , na.rm = TRUE)
MIS_complete_cell <- sum(distransform_MIS$MeanIntensity > 1.5 , na.rm = TRUE)
MIS_total_cell <- sum(distransform_MIS$MeanIntensity >= 0, na.rm = TRUE)

# Calculate percentages
MIS_incomplete_cell_percentage <- (MIS_incomplete_cell / MIS_total_cell) * 100
MIS_complete_cell_percentage <- (MIS_complete_cell  / MIS_total_cell) * 100

# Print results
cat("Incomplete Cell Number: ", MIS_incomplete_cell, "\n")
cat("Complete Cell Number: ", MIS_complete_cell, "\n")
cat("Total Cell Number: ", MIS_total_cell, "\n")
cat("Incomplete Cell Percentage: ", MIS_incomplete_cell_percentage, "%\n")
cat("Complete Cell Percentage: ", MIS_complete_cell_percentage, "%\n")

# VGP 
# Calculate counts for each condition
VGP_incomplete_cell <- sum(distransform_VGP$MeanIntensity <= 1.5 , na.rm = TRUE)
VGP_complete_cell <- sum(distransform_VGP$MeanIntensity > 1.5 , na.rm = TRUE)
VGP_total_cell <- sum(distransform_VGP$MeanIntensity >= 0, na.rm = TRUE)

# Calculate percentages
VGP_incomplete_cell_percentage <- (VGP_incomplete_cell / VGP_total_cell) * 100
VGP_complete_cell_percentage <- (VGP_complete_cell  / VGP_total_cell) * 100

# Print results
cat("Incomplete Cell Number: ", VGP_incomplete_cell, "\n")
cat("Complete Cell Number: ", VGP_complete_cell, "\n")
cat("Total Cell Number: ", VGP_total_cell, "\n")
cat("Incomplete Cell Percentage: ", VGP_incomplete_cell_percentage, "%\n")
cat("Complete Cell Percentage: ", VGP_complete_cell_percentage, "%\n")

