read_data <- function(input_files) {
  t_data <- Read10X(input_files)
  t_data <- CreateSeuratObject(
    counts = t_data,
    project = input_files,
    min.cells = 3
  )
  
  t_data$log10GenesPerUMI <-
    log10(t_data$nFeature_RNA) / log10(t_data$nCount_RNA)
  t_data$percent_mt <- PercentageFeatureSet(t_data, pattern = "^mt-")
  t_data$percent_ribo <- PercentageFeatureSet(t_data, pattern = "^Rps")
  
  return(t_data)
}

read_data_h5 <- function(input_files) {
  t_data <- Read10X_h5(input_files)
  t_data <- CreateSeuratObject(
    counts = t_data,
    project = input_files,
    min.cells = 3,
    min.features = 200
  )
  
  t_data$log10GenesPerUMI <-
    log10(t_data$nFeature_RNA) / log10(t_data$nCount_RNA)
  t_data$percent_mt <- PercentageFeatureSet(t_data, pattern = "^mt-")
  t_data$percent_ribo <- PercentageFeatureSet(t_data, pattern = "^Rps")
  
  return(t_data)
}

# nCount_RNA
n_count_rnaf <- function(t_combined) {
  t_combined %>%
    ggplot(aes(
      color = .$orig.ident, # nolint
      x = .$nCount_RNA,
      y = .$orig.ident,
      fill = .$orig.ident
    )) +
    geom_density_ridges(alpha = 0.2, scale = 10) +
    scale_x_log10() +
    theme_classic() +
    ylab("Cell density") +
    xlab("nCount_RNA (UMI)") +
    ggtitle("Normalized Count RNA") +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_vline(xintercept = c(2000, 31000))
}

# nFeature_RNA
n_feature_rnaf <- function(t_combined) {
  t_combined %>%
    ggplot(aes(
      color = .$orig.ident, # nolint
      x = .$nFeature_RNA,
      y = .$orig.ident,
      fill = .$orig.ident
    )) +
    geom_density_ridges(alpha = 0.2, scale = 10) +
    scale_x_log10() +
    theme_classic() +
    ylab("Cell density") +
    xlab("nFeature_RNA") +
    ggtitle("Normalized Feature RNA") +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_vline(xintercept = c(900, 6500))
}


# percent.mito
percent_mitof <- function(t_combined) {
  t_combined %>%
    ggplot(aes(
      color = .$orig.ident, # nolint
      x = .$percent_mt + 0.2,
      y = .$orig.ident,
      fill = .$orig.ident
    )) +
    geom_density_ridges(alpha = 0.2, scale = 10) +
    scale_x_log10() +
    theme_classic() +
    ylab("Cell density") +
    xlab("percent_mito + 0.2") +
    ggtitle("Mitochondrial RNA") +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_vline(xintercept = c(0.9, 12.5))
}


# percent.ribo
percent_ribof <- function(t_combined) {
  t_combined %>%
    ggplot(aes(
      color = .$orig.ident, # nolint
      x = .$percent_ribo,
      y = .$orig.ident,
      fill = .$orig.ident
    )) +
    geom_density_ridges(alpha = 0.2, scale = 10) +
    theme_classic() +
    ylab("Cell density") +
    xlab("percent_ribo") +
    ggtitle("Ribosomal RNA") +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_vline(xintercept = c(18))
}


# complexity equals (log10Genes/UMI)
complexityf <- function(t_combined) {
  t_combined %>%
    ggplot(aes(
      x = .$log10GenesPerUMI, # nolint
      color = .$orig.ident,
      y = .$orig.ident,
      fill = .$orig.ident
    )) +
    geom_density_ridges(alpha = 0.2, scale = 10) +
    theme_classic() +
    ylab("Cell density") +
    xlab("log10(nFeature_RNA/nCount_RNA)") +
    ggtitle("Complexity") +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_vline(xintercept = c(0.815, 0.935))
}

qc_functions <- function(t_combined) {
  t_combined %>%
    ggplot(aes(x = nCount_RNA, y = nFeature_RNA, color = percent_mt)) + # nolint
    geom_point() +
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method = lm) +
    scale_x_log10() +
    scale_y_log10() +
    theme_minimal() +
    geom_vline(xintercept = 500) +
    geom_hline(yintercept = 250) +
    facet_wrap(~orig.ident)
}


# REMOVE OUTLIERS
quality_control <- function(input_data) {
  subset <- input_data %>%
    subset(subset = nCount_RNA > 2000 & nCount_RNA < 31000 & # nolint
             nFeature_RNA > 900 & nFeature_RNA < 6500 & # nolint
             percent_mt > 0.9 & percent_mt < 12.5 & # nolint
             percent_ribo < 18 & # nolint
             log10GenesPerUMI > 0.815 & log10GenesPerUMI < 0.935) # nolint
  
  return(subset)
}

# Output a logical vector for every gene on whether
# the more than zero counts per cell
gene_filteringf <- function(filtered_seurat) {
  # Extract counts
  counts <- GetAssayData(object = filtered_seurat, slot = "counts")
  
  # Output a logical vector for every gene on whether
  # the more than zero counts per cell
  nonzero <- counts > 0
  
  # Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
  keep_genes <- Matrix::rowSums(nonzero) >= 10
  
  # Only keeping those genes expressed in more than 10 cells
  filtered_counts <- counts[keep_genes, ]
  
  # Reassign to filtered Seurat object
  filtered_seurat <- CreateSeuratObject(filtered_counts,
                                        meta.data = filtered_seurat@meta.data
  )
  return(NormalizeData(filtered_seurat))
}



# PRINT QC PLOT
view_qc_plots <- function(input_data) {
  plot1 <- FeatureScatter(input_data,
                          feature1 = "nCount_RNA", feature2 = "percent.mt"
  )
  plot2 <- FeatureScatter(input_data,
                          feature1 = "nCount_RNA", feature2 = "nFeature_RNA"
  )
  plot3 <- FeatureScatter(input_data,
                          feature1 = "nFeature_RNA", feature2 = "percent.mt"
  )
  plot4 <- FeatureScatter(input_data,
                          feature1 = "nFeature_RNA", feature2 = "percent.ribo"
  )
  plot1 + plot2 + plot3 + plot4
}


integration_pcaf <- function(a, b) {
  a_anchors <- FindIntegrationAnchors(
    object.list = list(a, b),
    dims = 1:50
  )
  a_combined <- IntegrateData(anchorset = a_anchors, dims = 1:50)
  
  DefaultAssay(a_combined) <- "integrated"
  a_combined <- ScaleData(a_combined)
  
  n_of_npcs <- 50
  a_combined <- RunPCA(a_combined, npcs = n_of_npcs)
  DimHeatmap(a_combined, dims = 1:5, ncol = 5) %>% print()
  ElbowPlot(a_combined, ndims = n_of_npcs) %>% print()
  
  ndims_number <- 30
  a_combined <- a_combined %>%
    FindNeighbors(reduction = "pca", dims = 1:ndims_number) %>%
    FindClusters(resolution = 0.8) %>%
    RunUMAP(reduction = "pca", dims = 1:ndims_number) %>%
    RunTSNE(dims = 1:ndims_number)
  
  print(DimPlot(a_combined, label = TRUE))
  Idents(a_combined) <- "orig.ident"
  Idents(a_combined) <- "seurat_clusters"
  return(a_combined)
}







