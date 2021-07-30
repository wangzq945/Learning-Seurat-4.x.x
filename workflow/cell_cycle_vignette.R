# Cell-Cycle Scoring and Regression ----

# package
library(Seurat)

# path
c("output") %>% walk(function(x){if(!dir.exists(x)) {dir.create(x)}})

## Load data ----

# Read in the expression matrix The first row is a header row, the first column is rownames
exp.mat <- read.table(file ="data/nestorawa_forcellcycle_expressionMatrix.txt", header = TRUE, as.is = TRUE, row.names = 1)

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

## Standard workflow ----

# Create our Seurat object and complete the initalization steps
marrow <- CreateSeuratObject(counts = exp.mat)
marrow <- NormalizeData(marrow)
marrow <- FindVariableFeatures(marrow, selection.method = "vst")
marrow <- ScaleData(marrow, features = rownames(marrow))

marrow <- RunPCA(marrow, features = VariableFeatures(marrow), ndims.print = 6:10, nfeatures.print = 10)
DimHeatmap(marrow, dims = c(8, 10))

## Assign Cell-Cycle Scores ----

marrow <- CellCycleScoring(marrow, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(marrow[[]])

# Visualize the distribution of cell cycle markers across
RidgePlot(marrow, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by phase
marrow <- RunPCA(marrow, features = c(s.genes, g2m.genes))
DimPlot(marrow)

## Regress out cell cycle scores during data scaling ----

marrow <- ScaleData(marrow, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(marrow))

# Now, a PCA on the variable genes no longer returns components associated with cell cycle
marrow <- RunPCA(marrow, features = VariableFeatures(marrow), nfeatures.print = 10)

# When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
marrow <- RunPCA(marrow, features = c(s.genes, g2m.genes))
DimPlot(marrow)

## Alternate Workflow ----

marrow$CC.Difference <- marrow$S.Score - marrow$G2M.Score
marrow <- ScaleData(marrow, vars.to.regress = "CC.Difference", features = rownames(marrow))

# cell cycle effects strongly mitigated in PCA
marrow <- RunPCA(marrow, features = VariableFeatures(marrow), nfeatures.print = 10)

# when running a PCA on cell cycle genes, actively proliferating cells remain distinct from G1 cells however, within actively proliferating cells, G2M and S phase cells group together
marrow <- RunPCA(marrow, features = c(s.genes, g2m.genes))
DimPlot(marrow)

## Save ----

saveRDS(marrow, file = "output/marrow.rds")
