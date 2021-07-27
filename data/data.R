# The input data ----

# package
library(Seurat)
library(SeuratData)

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
