# Using sctransform in Seurat ----

# package
library(Seurat)
library(ggplot2)
library(sctransform)

# path
c("output") %>% walk(function(x){if(!dir.exists(x)) {dir.create(x)}})

## SCT ----

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = paste0(mainPath, "/data/filtered_gene_bc_matrices/hg19/"))

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# store mitochondrial percentage in object meta data
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = 'percent.mt')

# run sctransform
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)

## Standard steps in the Seurat workflow for visualization and clustering ----

# Visualize canonical marker genes as violin plots.
VlnPlot(pbmc, features = c("CD8A", "GZMK", "CCL5", "S100A4", "ANXA1", "CCR7", "ISG15", "CD3D"), pt.size = 0.2, ncol = 4)

# Visualize canonical marker genes on the sctransform embedding.
FeaturePlot(pbmc, features = c("CD8A", "GZMK", "CCL5", "S100A4", "ANXA1", "CCR7"), pt.size = 0.2, ncol = 3)
FeaturePlot(pbmc, features = c("CD3D", "ISG15", "TCL1A", "FCER2", "XCL1", "FCGR3A"), pt.size = 0.2, ncol = 3)

## Save ----

saveRDS(pbmc, file = "output/pbmc3k_final_sct.rds")
