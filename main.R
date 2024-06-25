# LOAD PACKAGES
library(Seurat)
library(dplyr)
# library(tibble)
# library(cowplot)
library(ggplot2)
library(ggridges)
# library(EnhancedVolcano)
# library(gridExtra)
# library(SingleCellExperiment)
# library(pheatmap)
# library(patchwork)
# library(ggpubr)

setwd("C:/Users/heung/OneDrive/Documents/scData/15. gdTcell_PWH/")
source('functions.R')


# Preprocessing ================================================================
data1 <- read_data("PJH")
data2 <- read_data_h5("KCW.h5")

t_combined <- 
  rbind(data1@meta.data, data2@meta.data)


t_combined %>% n_count_rnaf()
t_combined %>% n_feature_rnaf()
t_combined %>% percent_mitof()
t_combined %>% percent_ribof()
t_combined %>% complexityf()

data1 <- data1 %>% quality_control() %>% NormalizeData()
data2 <- data2 %>% quality_control() %>% NormalizeData()

combined <- integration_pcaf(data1, data2)
saveRDS(combined, "combined.rds")

combined <- readRDS("combined.rds")

DefaultAssay(combined) <- "RNA"
DefaultAssay(combined) <- "Integrated"
VlnPlot(combined, "Trdc")


# Histogram 그리기 =============================================================
df <- FetchData(
  object = combined,
  vars = c("Ptprc", "Trac", "Trdc"),
)

# histogram
ggplot(df, aes(x = Ptprc))+ geom_histogram() + geom_vline(xintercept = 0.4) + theme_minimal()
ggplot(df, aes(x = Trac))+ geom_histogram() + geom_vline(xintercept = 0.3) + theme_minimal()
ggplot(df, aes(x = Trdc))+ geom_histogram() + geom_vline(xintercept = 0.2) + theme_minimal()

# ridge plot
ggplot(df, aes(x = Ptprc, y = stat(density))) +
  geom_density_ridges() + theme_minimal() + geom_vline(xintercept = 0.4)
ggplot(df, aes(x = Trac, y = stat(density))) +
  geom_density_ridges() + theme_minimal() + geom_vline(xintercept = 0.3)
ggplot(df, aes(x = Trdc, y = stat(density))) +
  geom_density_ridges() + theme_minimal() + geom_vline(xintercept = 0.2)


# gdT filtering ================================================================
Ptprc_high <- WhichCells(combined, expression = Ptprc > 0.4)
Trac_low <- WhichCells(combined, expression = Trac < 0.3)
Trdc_high <- WhichCells(combined, expression = Trdc > 0.2)

gdT <- intersect(intersect(Ptprc_high, Trac_low), Trdc_high)

# Seurat object cell metadata에 gdT column 추가 (False and True)
combined <- AddMetaData(
  combined,
  metadata = rep("False", length(Cells(combined))),
  col.name = "gdT"
)
combined@meta.data[gdT, "gdT"] <- "True"

# UMAP 확인인
DimPlot(combined, group.by = "gdT") +
  scale_color_manual(values = c("grey", "red"), guide = FALSE)


# gdT 만 뽑아서 gdT seurat object로 저장장 =====================================
gdT <- combined[, (combined$gdT %in% "True")]

VlnPlot(gdT, c("Ptprc", "Trac", "Trdc"), group.by = 'gdT')
