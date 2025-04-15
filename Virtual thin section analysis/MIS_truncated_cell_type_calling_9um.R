library(dplyr)

#MHCI
MHCI_9um <- read.csv("Z:/F8ll_results/In_situ_9_micron/MHCI_73_105_9um_truncated_stats.csv")
MHCI_36um <- read.csv("Z:/F8ll_results/F8iia-quantification5.csv")
Merge <- left_join(x = MHCI_9um, y = MHCI_36um, by = "CELLID")
write.csv(Merge, file = "Z:/F8ll_results/In_situ_9_micron/MHCI_73_105_9um_stats_volume_celltype.csv")

#LAG3
MHCI_9um <- read.csv("Z:/F8ll_results/In_situ_9_micron/MHCI_73_105_9um_truncated_stats.csv")
LAG3_9um <- read.csv("Z:/F8ll_results/In_situ_9_micron/Dataset1-LSP13626-invasive_margin_LAG3_73_105_9um_stats.csv")
Merge <- left_join(x = MHCI_9um, y = LAG3_9um, by = "CELLID")
write.csv(Merge, file = "Z:/F8ll_results/In_situ_9_micron/MHCI_73_105_9um_stats_LAG3_truncated_celltype.csv")

#MX1
MHCI_9um <- read.csv("Z:/F8ll_results/In_situ_9_micron/MHCI_73_105_9um_truncated_stats.csv")
MX1_9um <- read.csv("Z:/F8ll_results/In_situ_9_micron/Dataset1-LSP13626-melanoma_in-situ_MX1_73_105_9um_stats.csv")
Merge <- left_join(x = MHCI_9um, y = MX1_9um, by = "CELLID")
write.csv(Merge, file = "Z:/F8ll_results/In_situ_9_micron/MHCI_73_105_9um_stats_MX1_truncated_celltype.csv")

#GZMB
MHCI_9um <- read.csv("Z:/F8ll_results/In_situ_9_micron/MHCI_73_105_9um_truncated_stats.csv")
GZMB_9um <- read.csv("Z:/F8ll_results/In_situ_9_micron/Dataset1-LSP13626-invasive_margin_GZMB_73_105_9um_stats.csv")
Merge <- left_join(x = MHCI_9um, y = GZMB_9um, by = "CELLID")
write.csv(Merge, file = "Z:/F8ll_results/In_situ_9_micron/MHCI_73_105_9um_stats_GZMB_truncated_celltype.csv")

#CD103
MHCI_9um <- read.csv("Z:/F8ll_results/In_situ_9_micron/MHCI_73_105_9um_truncated_stats.csv")
CD103_9um <- read.csv("Z:/F8ll_results/In_situ_9_micron/Dataset1-LSP13626-melanoma_in-situ_CD103_73_105_9um_stats.csv")
Merge <- left_join(x = MHCI_9um, y = CD103_9um, by = "CELLID")
write.csv(Merge, file = "Z:/F8ll_results/In_situ_9_micron/MHCI_73_105_9um_stats_CD103_truncated_celltype.csv")

#MART1
MHCI_9um <- read.csv("Z:/F8ll_results/In_situ_9_micron/MHCI_73_105_9um_truncated_stats.csv")
MART1_9um <- read.csv("Z:/F8ll_results/In_situ_9_micron/Dataset1-LSP13626-melanoma_in-situ_MART1_73_105_9um_stats.csv")
Merge <- left_join(x = MHCI_9um, y = MART1_9um, by = "CELLID")
write.csv(Merge, file = "Z:/F8ll_results/In_situ_9_micron/MHCI_73_105_9um_stats_MART1_truncated_celltype.csv")

#PDL1
MHCI_9um <- read.csv("Z:/F8ll_results/In_situ_9_micron/MHCI_73_105_9um_truncated_stats.csv")
PDL1_9um <- read.csv("Z:/F8ll_results/In_situ_9_micron/Dataset1-LSP13626-melanoma_in-situ_PDL1_73_105_9um_stats.csv")
Merge <- left_join(x = MHCI_9um, y = PDL1_9um, by = "CELLID")
write.csv(Merge, file = "Z:/F8ll_results/In_situ_9_micron/MHCI_73_105_9um_stats_PDL1_truncated_celltype.csv")

#PD1
MHCI_9um <- read.csv("Z:/F8ll_results/In_situ_9_micron/MHCI_73_105_9um_truncated_stats.csv")
PD1_9um <- read.csv("Z:/F8ll_results/In_situ_9_micron/Dataset1-LSP13626-melanoma_in-situ_PD1_73_105_9um_stats.csv")
Merge <- left_join(x = MHCI_9um, y = PD1_9um, by = "CELLID")
write.csv(Merge, file = "Z:/F8ll_results/In_situ_9_micron/MHCI_73_105_9um_stats_PD1_truncated_celltype.csv")

#CD4
MHCI_9um <- read.csv("Z:/F8ll_results/In_situ_9_micron/MHCI_73_105_9um_truncated_stats.csv")
CD4_9um <- read.csv("Z:/F8ll_results/In_situ_9_micron/Dataset1-LSP13626-melanoma_in-situ_CD4_73_105_9um_stats.csv")
Merge <- left_join(x = MHCI_9um, y = CD4_9um, by = "CELLID")
write.csv(Merge, file = "Z:/F8ll_results/In_situ_9_micron/MHCI_73_105_9um_stats_CD4_truncated_celltype.csv")

#CD8
MHCI_9um <- read.csv("Z:/F8ll_results/In_situ_9_micron/MHCI_73_105_9um_truncated_stats.csv")
CD8_9um <- read.csv("Z:/F8ll_results/In_situ_9_micron/Dataset1-LSP13626-melanoma_in-situ_CD8_73_105_9um_stats.csv")
Merge <- left_join(x = MHCI_9um, y = CD8_9um, by = "CELLID")
write.csv(Merge, file = "Z:/F8ll_results/In_situ_9_micron/MHCI_73_105_9um_stats_CD8_truncated_celltype.csv")

#SOX10
MHCI_9um <- read.csv("Z:/F8ll_results/In_situ_9_micron/MHCI_73_105_9um_truncated_stats.csv")
SOX10_9um <- read.csv("Z:/F8ll_results/In_situ_9_micron/Dataset1-LSP13626-melanoma_in-situ_SOX10_73_105_9um_stats.csv")
Merge <- left_join(x = MHCI_9um, y = SOX10_9um, by = "CELLID")
write.csv(Merge, file = "Z:/F8ll_results/In_situ_9_micron/MHCI_73_105_9um_stats_SOX10_truncated_celltype.csv")

# CD8 PD1 double positive cells
CD8_PD1 <- read.csv("Z:/F8ll_results/In_situ_9_micron/MHCI_73_105_9um_truncated_stats_CD8_PD1_positive.csv")
CD8_9um <- read.csv("Z:/F8ll_results/In_situ_9_micron/Dataset1-LSP13626-melanoma_in-situ_CD8_73_105_9um_stats.csv")
Merge <- left_join(x = CD8_PD1, y = CD8_9um, by = "CELLID")
write.csv(Merge, file = "Z:/F8ll_results/In_situ_9_micron/MHCI_73_105_9um_stats_CD8_PD1_truncated_celltype_CD8.csv")

CD8_PD1 <- read.csv("Z:/F8ll_results/In_situ_9_micron/MHCI_73_105_9um_truncated_stats_CD8_PD1_positive.csv")
PD1_9um <- read.csv("Z:/F8ll_results/In_situ_9_micron/Dataset1-LSP13626-melanoma_in-situ_PD1_73_105_9um_stats.csv")
Merge <- left_join(x = CD8_PD1, y = PD1_9um, by = "CELLID")
write.csv(Merge, file = "Z:/F8ll_results/In_situ_9_micron/MHCI_73_105_9um_stats_CD8_PD1_truncated_celltype_PD1.csv")
