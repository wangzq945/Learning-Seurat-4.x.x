# Mapping and annotating query datasets ----

# package
library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)

# path
c("output") %>% walk(function(x){if(!dir.exists(x)) {dir.create(x)}})

## Load data ----

# human pancreatic islet cell datasets produced across four technologies, CelSeq (GSE81076) CelSeq2 (GSE85241), Fluidigm C1 (GSE86469), and SMART-Seq2 (E-MTAB-5061)
panc8 <- readRDS("data/panc8.rds")

## Dataset pre-processing ----

# split the combined object into a list, with each dataset as an element
pancreas.list <- SplitObject(panc8, split.by = "tech")
pancreas.list <- pancreas.list[c("celseq", "celseq2", "fluidigmc1", "smartseq2")]

# perform standard preprocessing (log-normalization), and identify variable features individually for each
for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

## Integration of 3 pancreatic islet cell datasets ----

# Assembly of multiple distinct scRNA-seq datasets into an integrated reference

# identify integration anchors 
reference.list <- pancreas.list[c("celseq", "celseq2", "smartseq2")]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)

# pass these anchors to integrate datasets
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)

# use this new integrated matrix for downstream analysis and visualization
# switch to integrated assay. The variable features of this assay are automatically
# set during IntegrateData
DefaultAssay(pancreas.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend()
p1 + p2

## Cell type classification using an integrated reference ----

# Transfer of cell type labels from a reference dataset onto a new query dataset
#   identify transfer anchors
#   pass these anchors to transfer label
#   add the predictions to the query metadata
pancreas.query <- pancreas.list[["fluidigmc1"]]
pancreas.anchors <- FindTransferAnchors(reference = pancreas.integrated, query = pancreas.query, dims = 1:30)
predictions <- TransferData(anchorset = pancreas.anchors, refdata = pancreas.integrated$celltype, dims = 1:30)
pancreas.query <- AddMetaData(pancreas.query, metadata = predictions)

# evaluate the prediction
pancreas.query$prediction.match <- pancreas.query$predicted.id == pancreas.query$celltype
table(pancreas.query$prediction.match)

# verify the prediction
table(pancreas.query$predicted.id)
VlnPlot(pancreas.query, c("REG1A", "PPY", "SST", "GHRL", "VWF", "SOX10"), group.by = "predicted.id")

## Unimodal UMAP Projection ----

pancreas.integrated <- RunUMAP(pancreas.integrated, dims = 1:30, reduction = "pca", return.model = TRUE)
pancreas.query <- MapQuery(
  anchorset = pancreas.anchors, 
  reference = pancreas.integrated,
  query = pancreas.query,
  refdata = list(celltype = 'celltype'),
  reference.reduction = 'pca',
  reduction.model = 'umap'
)

p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(pancreas.query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2
