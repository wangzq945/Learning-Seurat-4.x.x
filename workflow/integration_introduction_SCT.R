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

# SCT
ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform) 

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)

## Perform integration ----

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = 'SCT', anchor.features = features)

# this command creates an 'integrated' data assay
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = 'SCT')

## Perform an integrated analysis ----

# Now we can run a single integrated analysis on all cells!

# specify that we will perform downstream analysis on the corrected data
# note that the original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined.sct) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined.sct <- ScaleData(immune.combined.sct, verbose = FALSE)
immune.combined.sct <- RunPCA(immune.combined.sct, npcs = 30, verbose = FALSE)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:30)
immune.combined.sct <- FindNeighbors(immune.combined.sct, reduction = "pca", dims = 1:30)
immune.combined.sct <- FindClusters(immune.combined.sct, resolution = 0.5)

# Visualization
p1 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined.sct, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

## Save ----

saveRDS(immune.combined.sct, file = "output/immune.combined.sct.rds")
