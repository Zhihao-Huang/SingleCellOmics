.libPaths('C:\\Users\\Administrator\\Documents\\Rlibs\\Seurat3.1')
library(dplyr)
library(Seurat)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Load the PBMC dataset
# data was downloaded from https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
pbmc.data <- Read10X(data.dir = "./filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst",nfeatures = 2000)
#pbmc <- FindVariableFeatures(pbmc, selection.method = "mean.var.plot",mean.cutoff=c(0.0125,2),dispersion.cutoff=c(0.5,Inf))
length(VariableFeatures(pbmc))
pbmc@assays$RNA@var.features <- c(pbmc@assays$RNA@var.features,'CCR7')
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
# plot variable features with and without labels
#plot1 <- VariableFeaturePlot(pbmc,selection.method='mean.var.plot')
plot1 <- VariableFeaturePlot(pbmc)
LabelPoints(plot = plot1, points = 'CCR7', repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc <- RunTSNE(pbmc, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")
DimPlot(pbmc, reduction = "tsne")
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
                               "CD8A"))
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, ident.2 = 1, logfc.threshold = 0.25, test.use = "wilcox",only.pos = T)
FeaturePlot(pbmc, features = rownames(cluster0.markers)[10:17],ncol=3)

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5,repel = T) + NoLegend()
FeaturePlot(pbmc, features = c('CCR7','CD4','CD8A','S100A4'),pt.size = 0.5)


###select PCs in which SELL and S100A4 are more variable
pcmat <- pbmc@reductions$pca@feature.loadings
allgene_order <- apply(abs(pcmat),2,order)
gene_pos <- which(rownames(pcmat) %in% c('CCR7','S100A4'))
gene_order <- allgene_order[gene_pos,]
select_PCs <- colnames(pcmat)[apply(gene_order,2,function(x){sum(x) < 700})]
select_PCs <- as.numeric(gsub('PC_','',select_PCs))
select_PCs
#[1]  2  5  6  7 26 30 32
pcmat <- pbmc@reductions$pca@feature.loadings
allgene_order <- apply(abs(pcmat),2,order)
gene_pos <- which(rownames(pcmat) %in% c('PPBP'))
gene_order2 <- allgene_order[gene_pos,]
select_PCs2 <- colnames(pcmat)[gene_order2 < 700]
select_PCs2 <- as.numeric(gsub('PC_','',select_PCs2))
select_PCs2 <- select_PCs2[3< select_PCs2 & select_PCs2 < 10]
select_PCs2
#[1] 4 8
pbmc <- RunUMAP(pbmc, dims = unique(c(select_PCs,select_PCs2)))
DimPlot(pbmc, reduction = "umap",pt.size = 0.5)

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
FeaturePlot(pbmc, features = c('CCR7','CD4','CD8A','S100A4'),pt.size = 0.5)


#####SVD test
svdresult <- svd(t(pbmc@assays$RNA@scale.data[VariableFeatures(pbmc),]))
feature.loadings <- svdresult$V
cell.embeddings <- svdresult$u %*% diag(svdresult$d)
