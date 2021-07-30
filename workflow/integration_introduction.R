# Introduction to scRNA-seq integration ----

# package
library(Seurat)
library(patchwork)

# path
c("output") %>% walk(function(x){if(!dir.exists(x)) {dir.create(x)}})

## Load data ----

ifnb <- readRDS(paste0(mainPath, "/data/ifnb.rds"))

## Setup the Seurat objects ----

# split the dataset into a list of two seurat objects (stim and CTRL)
ifnb.list <- SplitObject(ifnb, split.by = "stim")

# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list)

## Perform integration ----

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)

# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)

## Perform an integrated analysis ----

# Now we can run a single integrated analysis on all cells!

# specify that we will perform downstream analysis on the corrected data
# note that the original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

## Save ----

saveRDS(immune.combined, file = "output/immune.combined.rds")
