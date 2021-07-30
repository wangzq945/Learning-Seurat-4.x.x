# The input data ----

# package
library(Seurat)
library(SeuratData)

# for: integration_mapping
# human pancreatic islet cell datasets produced across four technologies, CelSeq (GSE81076) CelSeq2 (GSE85241), Fluidigm C1 (GSE86469), and SMART-Seq2 (E-MTAB-5061)
InstallData("panc8")
data("panc8")
saveRDS(panc8, "data/panc8.rds")

# for: spatial_vignette
# dataset of sagital mouse brain slices generated using the Visium v1 chemistry
# There are two serial anterior sections, and two (matched) serial posterior sections.
InstallData("stxBrain")
brain <- LoadData("stxBrain", type = "anterior1")
saveRDS(brain, "data/stxBrain.anterior1.rds")
brain2 <- LoadData("stxBrain", type = "posterior1")
saveRDS(brain2, "data/stxBrain.posterior1.rds")

# for: spatial_vignette
# dataset generated using Slide-seq v2 of the mouse hippocampus
InstallData("ssHippo")
slide.seq <- LoadData("ssHippo")
saveRDS(slide.seq, "data/ssHippo.rds")
