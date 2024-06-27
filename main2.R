
setwd("C:/")

# Read data=====================================================================
data1 <- Read10X("PJH")
data1 <- CreateSeuratObject(counts = data1, project = "WT1", min.cells = 3, min.features = 200)

data2 <- Read10X_h5("KCW.h5")
data2 <- CreateSeuratObject(counts = data2, project = "WT2", min.cells = 3, min.features = 200)

# Quality Check=================================================================
data1$pct_mt <- PercentageFeatureSet(data1, pattern = "^mt-")
data2$pct_mt <- PercentageFeatureSet(data2, pattern = "^mt-")

VlnPlot(data1, features = c("nFeature_RNA", "nCount_RNA", "pct_mt"))
VlnPlot(data2, features = c("nFeature_RNA", "nCount_RNA", "pct_mt"))

plot1 <- FeatureScatter(data1, feature1 = "nCount_RNA", feature2 = "pct_mt")
plot2 <- FeatureScatter(data1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(data1, feature1 = "nFeature_RNA", feature2 = "pct_mt")
plot1 + plot2 + plot3

plot1 <- FeatureScatter(data2, feature1 = "nCount_RNA", feature2 = "pct_mt")
plot2 <- FeatureScatter(data2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(data2, feature1 = "nFeature_RNA", feature2 = "pct_mt")
plot1 + plot2 + plot3

# Quality Control===============================================================
data1 <- subset(data1, subset = pct_mt < 10 & nCount_RNA < 40000 & nFeature_RNA > 400)
data2 <- subset(data2, subset = pct_mt < 10 & nCount_RNA < 40000 & nFeature_RNA > 400)

# Quality Control 됐는지 확인 ==================================================
plot1 <- FeatureScatter(data1, feature1 = "nCount_RNA", feature2 = "pct_mt")
plot2 <- FeatureScatter(data1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(data1, feature1 = "nFeature_RNA", feature2 = "pct_mt")
plot1 + plot2 + plot3

plot1 <- FeatureScatter(data2, feature1 = "nCount_RNA", feature2 = "pct_mt")
plot2 <- FeatureScatter(data2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(data2, feature1 = "nFeature_RNA", feature2 = "pct_mt")
plot1 + plot2 + plot3

# Normalization ================================================================
data1 <- NormalizeData(data1)
data2 <- NormalizeData(data2)

# Clustering ===================================================================
data1 <- FindVariableFeatures(data1)
data1 <- ScaleData(data1)
data1 <- RunPCA(data1, npcs = 50)
DimHeatmap(data1, dims = 1:5, ncol = 5)
ElbowPlot(data1, ndims = 50)

data1 <- FindNeighbors(data1, reduction = "pca", dims = 30)
data1 <- FindClusters(data1, resolution = 1)
data1 <- RunUMAP(data1, reduction = "pca", dims = 1:30)
data1 <- RunTSNE(data1, dims = 1:30)
DimPlot(data1, reduction = "umap")
DimPlot(data1, reduction = "tsne")


data2 <- FindVariableFeatures(data2)
data2 <- ScaleData(data2)
data2 <- RunPCA(data2, npcs = 50)
DimHeatmap(data2, dims = 1:5, ncol = 5)
ElbowPlot(data2, ndims = 50)

data2 <- FindNeighbors(data2, reduction = "pca", dims = 30)
data2 <- FindClusters(data2, resolution = 1)
data2 <- RunUMAP(data2, reduction = "pca", dims = 1:30)
data2 <- RunTSNE(data2, dims = 1:30)
DimPlot(data2, reduction = "umap")
DimPlot(data2, reduction = "tsne")

# Integration ==================================================================
a.anchors <- FindIntegrationAnchors(object.list = list(data1, data2), dims = 1:50)
a.combined <- IntegrateData(anchorset = a.anchors, dims = 1:50)

# Clustering the integrated data ===============================================
DefaultAssay(a.combined) <- "integrated"

a.combined <- ScaleData(a.combined)
a.combined <- RunPCA(a.combined, npcs = 50)
ElbowPlot(a.combined)
a.combined <- FindNeighbors(a.combined, reduction = 'pca', dims = 1:30)
a.combined <- FindClusters(a.combined, resolution = 0.5)
a.combined <- RunUMAP(a.combined, reduction = "pca", dims = 1:30)
a.combined <- RunTSNE(a.combined, dims = 1:40)

DimPlot(a.combined)

# Save and load the integrated file ============================================
saveRDS(a.combined, "a.combined.rds")

a.combined <- readRDS("a.combined.rds")


# ==============================================================================
# //                                                                          //
# //                       extract gd T cells                                 //
# //                                                                          //
# ==============================================================================

# Histogram 그리기==============================================================
DefaultAssay(a.combined) <- "RNA"

df <- FetchData(
  object = a.combined,
  vars = c("Ptprc", "Trac", "Trdc")
)

# histogram
ggplot(df, aes(x = Ptprc))+ geom_histogram() + geom_vline(xintercept = 0.3) + theme_minimal()
ggplot(df, aes(x = Trac))+ geom_histogram() + geom_vline(xintercept = 0.2) + theme_minimal()
ggplot(df, aes(x = Trdc))+ geom_histogram() + geom_vline(xintercept = 0.2) + theme_minimal()

# ridge plot
ggplot(df, aes(x = Ptprc, y = stat(density))) +
  geom_density_ridges() + theme_minimal() + geom_vline(xintercept = 0.3)
ggplot(df, aes(x = Trac, y = stat(density))) +
  geom_density_ridges() + theme_minimal() + geom_vline(xintercept = 0.2)
ggplot(df, aes(x = Trdc, y = stat(density))) +
  geom_density_ridges() + theme_minimal() + geom_vline(xintercept = 0.2)


# gd T filtering ===============================================================
Ptprc_high <- WhichCells(a.combined, expression = Ptprc > 0.3)
Trac_low <- WhichCells(a.combined, expression = Trac < 0.2)
Trdc_high <- WhichCells(a.combined, expression = Trdc > 0.2)

gdT <- intersect(intersect(Ptprc_high, Trac_low), Trdc_high)

# Seurat object cell metadata에 gdT column 추가 (False and True)
a.combined <- AddMetaData(
  a.combined,
  metadata = rep("False", length(Cells(a.combined))),
  col.name = "gdT"
)
a.combined@meta.data[gdT, "gdT"] <- "True"

# UMAP 확인
DimPlot(a.combined, group.by = "gdT") +
  scale_color_manual(values = c("grey", "red"), guide = FALSE)

a.combined@meta.data

# gdT 만 뽑아서 gdT seurat object로 저장 =======================================
gdT <- a.combined[, (a.combined$gdT %in% "True")]

Idents(gdT) <- "gdT"
VlnPlot(gdT, c("Ptprc", "Trac", "Cd3e"))


'''
cytokine_receptors <- list(
  # Il-1 Family Receptors
  # not found : Il1rapl1, Il1rapl2
  il_1_family = c("Il18r1", "Il18rap", "Il1r1", "Il1r2", "Il1rl1", "Sigirr"),

  # TNF receptors
  # not found : tnfrsf13c, Eda2r, Edar
  tnf_family = c(
    "Tnfrsf9", "Cd27", "Cd40", "Fas", "Tnfrsf21", "Pglyrp1",
    "Relt", "Tnfrsf1a", "Tnfrsf1b", "Tnfrsf11a", "Tnfrsf11b",
    "Tnfrsf12a", "Tnfrsf13b", "Tnfrsf14", "Tnfrsf17", "Tnfrsf18",
    "Tnfrsf19", "Tnfrsf25", "Ltbr", "Tnfrsf4", "Tnfrsf8",
    "Tnfrsf10b"
  ),
  for_ppts = c(
    "Tnfrsf21", "Tnfrsf1a", "Tnfrsf1b", "Tnfrs11a", "Tnfrsf13b",
    "Ltbr", "Ifnar1", "Ifnar2", "Ifngr1", "Ifngr2", "Il6ra", "Lifr",
    "Il10ra", "Il10rb", "Tgfbr1", "Tgfbr2", "Csf1r", "Il17ra", "Il4ra"
  ),

  # Interferon receptors
  interferons = c("Ifnar1", "Ifnar2", "Ifngr1", "Ifngr2"),

  # IL-6 receptors
  il_6_family = c(
    "Cntfr", "Il11ra1", "Il6ra", "Lepr", "Lifr", "Osmr",
    "Il31ra"
  ),

  # IL-10 receptors
  # not found : Il20ra, Il22ra2, Il22ra1
  il_10_family = c("Il10ra", "Il10rb", "Il20rb"),

  # TGF beta receptors
  tgf_beta_family = c("Acvr1c", "Atf2", "Eng", "Tgfbr1", "Tgfbr2", "Tgfbr3"),

  # other hot cytokine receptors
  # not found : Epha3, Vegfr2
  others = c(
    "Acvr2a", "Acvr2b", "Csf1r", "Erbb2", "Erbb3", "Egfr", "Igf1r",
    "Il17ra", "Il2ra", "Il3ra", "Lilrb4a", "Insr", "Il4ra", "Lilra6",
    "Il7r", "Ldlr", "Il2rb", "Pdgfrb"
  )
)
'''









