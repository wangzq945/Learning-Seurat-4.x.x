# Data visualization methods in Seurat ----

# package
library(Seurat)
library(ggplot2)
library(patchwork)

# path
c("output") %>% walk(function(x){if(!dir.exists(x)) {dir.create(x)}})

## Load data ----

pbmc3k.final <- readRDS("data/pbmc3k.final.rds")
pbmc3k.final$groups <- sample(c('group1', 'group2'), size = ncol(pbmc3k.final), replace = TRUE)
features <- c("LYZ", "CCL5", "IL32", "PTPRCAP", "FCGR3A", "PF4")
pbmc3k.final

## Five visualizations of marker feature expression ----

# Ridge plots - from ggridges. Visualize single cell expression distributions in each cluster
RidgePlot(pbmc3k.final, features = features, ncol = 2)

# Violin plot - Visualize single cell expression distributions in each cluster
VlnPlot(pbmc3k.final, features = features)

# Feature plot - visualize feature expression in low-dimensional space
FeaturePlot(pbmc3k.final, features = features)

# Dot plots - the size of the dot corresponds to the percentage of cells expressing the feature in each cluster. The color represents the average expression level
DotPlot(pbmc3k.final, features = features) + RotatedAxis()

# Single cell heatmap of feature expression
DoHeatmap(subset(pbmc3k.final, downsample = 100), features = features, size = 3)

## New additions to FeaturePlot ----

# Plot a legend to map colors to expression levels
FeaturePlot(pbmc3k.final, features = "MS4A1")

# Adjust the contrast in the plot
FeaturePlot(pbmc3k.final, features = "MS4A1", min.cutoff = 1, max.cutoff = 3)

# Calculate feature-specific contrast levels based on quantiles of non-zero expression. Particularly useful when plotting multiple markers
FeaturePlot(pbmc3k.final, features = c("MS4A1", "PTPRCAP"), min.cutoff = "q10", max.cutoff = "q90")

# Visualize co-expression of two features simultaneously
FeaturePlot(pbmc3k.final, features = c("MS4A1", "CD79A"), blend = TRUE)

# Split visualization to view expression by groups (replaces FeatureHeatmap)
FeaturePlot(pbmc3k.final, features = c("MS4A1", "CD79A"), split.by = "groups")

## Updated and expanded visualization functions ----

# Violin plots can also be split on some variable. Simply add the splitting variable to object metadata and pass it to the split.by argument
VlnPlot(pbmc3k.final, features = "percent.mt", split.by = "groups")

# SplitDotPlotGG has been replaced with the `split.by` parameter for DotPlot
DotPlot(pbmc3k.final, features = features, split.by = "groups") + RotatedAxis()

# SplitDotPlotGG has been replaced with the `split.by` parameter for DotPlot
DotPlot(pbmc3k.final, features = features, split.by = "groups") + RotatedAxis()

# DoHeatmap now shows a grouping bar, splitting the heatmap into groups or clusters. This can be changed with the `group.by` parameter
DoHeatmap(pbmc3k.final, features = VariableFeatures(pbmc3k.final)[1:100], cells = 1:500, size = 4, angle = 90) + NoLegend()

## Applying themes to plots ----

baseplot <- DimPlot(pbmc3k.final, reduction = "umap")

# Add custom labels and titles
baseplot + labs(title = "Clustering of 2,700 pbmc3k.finals")

# Use community-created themes, overwriting the default Seurat-applied theme Install ggmin with remotes::install_github('sjessa/ggmin')
baseplot + ggmin::theme_powerpoint()

# Seurat also provides several built-in themes, such as DarkTheme; for more details see ?SeuratTheme
baseplot + DarkTheme()

# Chain themes together
baseplot + FontSize(x.title = 20, y.title = 20) + NoLegend()

## Interactive plotting features ----

# Include additional data to display alongside cell names by passing in a data frame of information Works well when using FetchData
plot <- FeaturePlot(pbmc3k.final, features = "MS4A1")
HoverLocator(plot = plot, information = FetchData(pbmc3k.final, vars = c("ident", "PC_1", "nFeature_RNA")))

# using CellSelector
pbmc3k.final <- RenameIdents(pbmc3k.final, DC = "CD14+ Mono")
plot <- DimPlot(pbmc3k.final, reduction = "umap")
select.cells <- CellSelector(plot = plot)
head(select.cells)

Idents(pbmc3k.final, cells = select.cells) <- "NewCells"

# Now, we find markers that are specific to the new cells, and find clear DC markers
newcells.markers <- FindMarkers(pbmc3k.final, ident.1 = "NewCells", ident.2 = "CD14+ Mono", min.diff.pct = 0.3, only.pos = TRUE)
head(newcells.markers)

## Plotting Accessories ----

# LabelClusters and LabelPoints will label clusters (a coloring variable) or individual points on a ggplot2-based scatter plot
plot <- DimPlot(pbmc3k.final, reduction = "pca") + NoLegend()
LabelClusters(plot = plot, id = "ident")

# Both functions support `repel`, which will intelligently stagger labels and draw connecting lines from the labels to the points or clusters
LabelPoints(plot = plot, points = TopCells(object = pbmc3k.final[["pca"]]), repel = TRUE)

plot1 <- DimPlot(pbmc3k.final)
plot2 <- FeatureScatter(pbmc3k.final, feature1 = "LYZ", feature2 = "CCL5")
# Combine two plots
plot1 + plot2
# Remove the legend from all plots
(plot1 + plot2) & NoLegend()
