##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)

#### Load Packages ####
if (!require('tidyverse')) {install.packages('tidyverse'); library(tidyverse)}
if (!require('zellkonverter')) {install.packages('zellkonverter', repos = "http://bioconductor.org/packages/release/bioc"); library(zellkonverter)}
if (!require('SingleCellExperiment')) {install.packages('SingleCellExperiment', repos = "http://bioconductor.org/packages/release/bioc"); library(SingleCellExperiment)}
if (!require('Seurat')) {install.packages('Seurat'); library(Seurat)}
if (!require('reticulate')) {install.packages('reticulate'); library(reticulate)}
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require('zellkonverter')) {BiocManager::install("zellkonverter")}; library(zellkonverter)

#### Load Data ####
file_path <- "C:/Charlene/Other/HIT/HIT A-20250111T231827Z-002/"
Set_FileName <- "Crunch2_scRNAseq"
sce <- zellkonverter::readH5AD(paste0(file_path, Set_FileName, ".h5ad"))

# Extract counts and metadata
counts <- assay(sce, "X")
metadata <- colData(sce)


# Create Seurat object
seurat_object <- CreateSeuratObject(counts = counts, meta.data = as.data.frame(metadata))

#### Data Preprocessing ####
seurat_object <- NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
seurat_object <- ScaleData(seurat_object)
seurat_object <- RunPCA(seurat_object)

# Elbow Plot
ElbowPlot(seurat_object)

# Clustering and UMAP
seurat_object <- FindNeighbors(seurat_object, dims = 1:30)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)
seurat_object <- RunUMAP(seurat_object, dims = 1:30)

# UMAP Plot
DimPlot(seurat_object, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(seurat_object, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "status") + NoLegend()

DimPlot(seurat_object, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "status") +
  DimPlot(seurat_object, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "study")

DimPlot(seurat_object, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "individual")
DimPlot(seurat_object, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "annotation") 

# 保存處理後的 Seurat 對象
saveRDS(seurat_obj, file = "C:/Charlene/Other/HIT/Processed_Crunch2_seurat_obj.rds")
