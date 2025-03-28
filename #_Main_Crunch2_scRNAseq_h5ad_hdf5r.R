# 安裝並載入所需套件
if(!require("Seurat")) {install.packages("Seurat"); library(Seurat)}
if(!require("hdf5r")) {install.packages("hdf5r"); library(hdf5r)}
if(!require("SeuratDisk")) {install.packages("SeuratDisk"); library(SeuratDisk)}

# 檢查 h5ad 文件的內容
file_path <- "C:/Charlene/Other/HIT/HIT A-20250111T231827Z-002/Crunch2_scRNAseq.h5ad"
h5_file <- H5File$new(file_path, mode = "r")
print(h5_file)
print(h5_file[["obs"]])  # 查看 obs 部分
print(h5_file[["var"]])  # 查看 var 部分
h5_file$close()

# 有時候 .h5ad 文件的元數據格式與 Seurat 期望的結構不完全一致。嘗試直接轉換而不指定 assay，以避免潛在的錯誤。
Convert(file_path, dest = "h5seurat", overwrite = TRUE)

# 讀取轉換後的文件
seurat_file <- sub(".h5ad$", ".h5seurat", file_path)
seurat_obj <- LoadH5Seurat(seurat_file)


# 確認元數據並修正類別型欄位
if(!all(c("annotation", "individual", "status", "dysplasia") %in% colnames(seurat_obj@meta.data))) {
  stop("Missing required metadata columns")
}

# 手動處理類別型資料（若必要）
seurat_obj@meta.data <- seurat_obj@meta.data %>%
  mutate(annotation = as.factor(annotation),
         individual = as.factor(individual),
         status = as.factor(status),
         dysplasia = as.factor(dysplasia))


# 初步前處理
# 過濾低質量細胞（例如，基於特定基因數量或線粒體比例）
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# 資料標準化
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# 鑑別高變基因
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# 數據縮放
seurat_obj <- ScaleData(seurat_obj)

# 主成分分析 (PCA)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

# 最近鄰網絡和分群分析
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# 非線性降維 (UMAP)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# 繪製 UMAP 圖檢視結果
DimPlot(seurat_obj, reduction = "umap", group.by = "annotation")

# 保存處理後的 Seurat 對象
saveRDS(seurat_obj, file = "C:/Charlene/Other/HIT/Processed_Crunch2_seurat_obj.rds")

# 完成
