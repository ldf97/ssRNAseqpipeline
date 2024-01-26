setwd("...\Data")
#Loading Libraries
library(Seurat)
library(tidyverse)
# Loading the NSCLC dataset
#NSCLC == Non Small Cell Lung Cancer

nsclc.sparse.m <- Read10X_h5(filename='20k_NSCLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5')
#Determining what modalities are present. It is Gene Expression, Antibody Capture, Multiplex Capture. Choosing Gene expression
str(nsclc.sparse.m)
cts <- nsclc.sparse.m$'Gene Expression'

#Initialize the seurat object with the raw (non-normalized data)
#We are keeping all features that are expressed in least 3 cells, keep all cells with at least 200 features/genes
nsclc.seurat.obj <- CreateSeuratObject(counts = cts, project = "NSCLC",min.cells = 3, min.features = 200)
# Inspect seurat object
str(nsclc.seurat.obj)
nsclc.seurat.obj

# Step 1: Quality Control. Filter out low quality cells
#Metrics for filtering cells: number of features/genes cells, number of total molecules (ncount)
#Poor quality cell: low number of genes. Extremely high number of genes: possibly due to doublets, multiple cells sequenced together
#This will hamper downstream processing, so filtering.
#Percentage of mitochondrial genes--> Dying cells have higher % of mitochondrial genes

View(nsclc.seurat.obj@meta.data)

nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(nsclc.seurat.obj, pattern="^MT-")
View(nsclc.seurat.obj@meta.data)

#viewing percentage of mitochondrial DNA as violin plot
VlnPlot(nsclc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
#Plot two metrics on two different axes, in this case Count of RNA and Count of Genes.
#Good quality cells should have linear relationship between number of genes and number of RNA molecules
FeatureScatter(nsclc.seurat.obj, feature1="nCount_RNA", feature2="nFeature_RNA") +
               geom_smooth(method='lm')
#From plot: Lower right quadrant means that experiment resulted in the same cells being sequenced multiple times, resulting in duplicated RNA counts
#The upper left quadrant represents a high number of unique genes counted, but a lack of depth--meaning there were not enough cells with those sets of genes so they are low quality for downstream analysis
# The next step is filtering out low quality cells

# 2. Filtering--------
nsclc.seurat.obj <- subset(nsclc.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
                             percent.mt < 5)
#Now, we have cells with a higher number of genes that are between 200 and 2500, as well as low number of mitochondrial DNA percentages
# We could also filter out Ribosomal Genes, this depends on the requirements and expected biological outcome
#DoubletFinder is a package that filters out doublets--two cells clumped together and sequenced together and labelled as one cell

# 3. Normalize data ----------
#In order to compare gene expression across multiple cells we must normalize data
#Divide gene expression measurements in each cell by total expresion, multiply by a scaling factor, then log transform it
nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj)
str(nsclc.seurat.obj) #We can look at commands that are run on seurat object

# 4. Identify Highly Variable Features---------
# We only want to select a few cells with high cell-to-cell variations
# This has been found that focusing on this subset of genes in downstream analysis highlights biological signal in single cell data set
# We can use a function called "Find Variable Features"
nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, selection.method = "vst", nfeatures=2000) #we can increase number of features if desired

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(nsclc.seurat.obj),10)
#Plot variable features with and without labels
plot1 <- VariableFeaturePlot(nsclc.seurat.obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)


#### Debug area:
# Load the NSCLC dataset
nsclc.sparse.m <- Read10X_h5(filename = '20k_NSCLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5')
str(nsclc.sparse.m)
cts <-  nsclc.sparse.m$`Gene Expression`



# Initialize the Seurat object with the raw (non-normalized data).
nsclc.seurat.obj <- CreateSeuratObject(counts = cts, project = "NSCLC", min.cells = 3, min.features = 200)
str(nsclc.seurat.obj)
nsclc.seurat.obj
# 29552 features across 42081 samples



# 1. QC -------
View(nsclc.seurat.obj@meta.data)
# % MT reads
nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(nsclc.seurat.obj, pattern = "^MT-")
View(nsclc.seurat.obj@meta.data)

VlnPlot(nsclc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

# 2. Filtering -----------------
nsclc.seurat.obj <- subset(nsclc.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
                             percent.mt < 5)

# 3. Normalize data ----------
nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)

str(nsclc.seurat.obj)
# 4. Identify highly variable features --------------
nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(nsclc.seurat.obj), 10)

dev.new()
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(nsclc.seurat.obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

# 5. Scale Data------------------
# Single cell dataset has many sources of variation.
# Biological sources; difference in cell cycle
# Technical sources: batch effects
# Account for these biases so data doesn't cluster incorrectly
all.genes <- rownames(nsclc.seurat.obj)
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj, features = all.genes)

str(nsclc.seurat.obj)
# 6. Perform Linear Dimensional Reduction--------------------
nsclc.seurat.obj <- RunPCA(nsclc.seurat.obj, features = VariableFeatures(object = nsclc.seurat.obj))

# Visualize PCA results
print(nsclc.seurat.obj[["pca"]], dims= 1:5, nfeatures = 5)
DimHeatmap(nsclc.seurat.obj, dims = 1, cells =500, balanced = TRUE)

#Determine the dimensionality of data---choosing only statistically significant components
ElbowPlot(nsclc.seurat.obj)
# we want to capture all of the components that explain a higher proportion of the variance
# variance explained stops changing around 10-15 components will be chosen in this data set

# 7. Clustering ----------
nsclc.seurat.obj <- FindNeighbors(nsclc.seurat.obj, dims = 1:15) #provide number of principle components that capture most variation in data set

# Assign cells to clusters
# Understanding resolution
nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj, resolution= c(0.1, 0.3, 0.5,0.7, 1))
View(nsclc.seurat.obj@meta.data)

DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.1", label=TRUE)

# Setting up identity of clusters
Idents(nsclc.seurat.obj)
Idents(nsclc.seurat.obj) <- "RNA_snn_res.0.1"
# all cells identified by eight clusters

# non-linear dimensionality reduction------------
nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj, dims=1:15)

DimPlot(nsclc.seurat.obj, reduction ="umap")
