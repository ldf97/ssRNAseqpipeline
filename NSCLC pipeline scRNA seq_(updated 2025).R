# Set working directory
setwd("C:/Users/Dante/OneDrive - Drexel University/Documents/Genetics")

# Load Libraries
library(Seurat)
library(tidyverse)

# 1. Load NSCLC Dataset -------------------------------
# NSCLC == Non Small Cell Lung Cancer
nsclc.sparse.m <- Read10X_h5("20k_NSCLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5")
cts <- nsclc.sparse.m$'Gene Expression'

# 2. Initialize Seurat Object -------------------------
nsclc.seurat.obj <- CreateSeuratObject(
  counts = cts,
  project = "NSCLC",
  min.cells = 3,
  min.features = 200
)

# 3. Quality Control -----------------------------------
# Add % mitochondrial gene expression
nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(nsclc.seurat.obj, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(nsclc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

# 4. Filter Cells --------------------------------------
nsclc.seurat.obj <- subset(nsclc.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# 5. Normalize Data ------------------------------------
nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)

# 6. Find Highly Variable Features ---------------------
nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, selection.method = "vst", nfeatures = 2000)

# Visualize top 10 highly variable genes
top10 <- head(VariableFeatures(nsclc.seurat.obj), 10)
plot1 <- VariableFeaturePlot(nsclc.seurat.obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

# 7. Scale Data ----------------------------------------
all.genes <- rownames(nsclc.seurat.obj)
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj, features = all.genes)

# 8. Run PCA -------------------------------------------
nsclc.seurat.obj <- RunPCA(nsclc.seurat.obj, features = VariableFeatures(object = nsclc.seurat.obj))

# Visualize PCA
print(nsclc.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(nsclc.seurat.obj, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(nsclc.seurat.obj)

# Re-run neighbors and clustering at one resolution
nsclc.seurat.obj <- FindNeighbors(nsclc.seurat.obj, dims = 1:15)
nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj, resolution = 0.5)

# Now you should see the identity column
colnames(nsclc.seurat.obj@meta.data)
table(Idents(nsclc.seurat.obj))

# Visualize clusters (e.g., resolution 0.5)
DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.5", label = TRUE)
Idents(nsclc.seurat.obj) <- "RNA_snn_res.0.5"
table(Idents(nsclc.seurat.obj))

# 10. Non-linear Dimensionality Reduction (UMAP) -------
nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj, dims = 1:15)
DimPlot(nsclc.seurat.obj, reduction = "umap")

# 11. Biological Interpretation -- Cluster Annotation
# Find marker genes for each cluster
cluster.markers <- FindAllMarkers(nsclc.seurat.obj,
                                  only.pos = TRUE,
                                  min.pct = 0.25,
                                  logfc.threshold = 0.25)

# Top 5 markers per cluster
cluster.markers %>% group_by(cluster) %>% slice_max(n = 5, order_by = avg_log2FC)

# Visualize markers
FeaturePlot(nsclc.seurat.obj, features = c("EPCAM", "CD3D", "CD14", "MS4A1", "MKI67"), min.cutoff = "q9")


new.cluster.ids <- c("Epithelial cells",   # 0
                     "T cells",            # 1
                     "Monocytes",          # 2
                     "B cells",            # 3
                     "Proliferating cells",# 4
                     "Unknown",            # 5
                     "Unknown",            # 6
                     "Unknown",            # 7
                     "Unknown",            # 8
                     "Unknown",            # 9
                     "Unknown",            #10
                     "Unknown",            #11
                     "Unknown")            #12

# Now this will work:
names(new.cluster.ids) <- levels(nsclc.seurat.obj)
nsclc.seurat.obj <- RenameIdents(nsclc.seurat.obj, new.cluster.ids)

#Visualizing clusters
DimPlot(nsclc.seurat.obj, label = TRUE, pt.size = 0.5) + NoLegend()

