library(Matrix)
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(viridis))
library(glmGamPoi)
library(RColorBrewer)
options(future.globals.maxSize = 10000 * 1024^2)
setwd("../../snrna/")

# Predict cell types from references

referencehippocampus = readRDS("ref/external_data/brain_atlas_subsampled_sct_pca.rds")

DefaultAssay(combined.sct) <- "SCT"
DefaultAssay(referencehippocampus) <- "SCT"

transfer_anchors <- FindTransferAnchors(
    reference = referencehippocampus,
    query = combined.sct,
    reference.assay = "SCT",
    normalization.method = "SCT",
    reference.reduction = "pca",
    recompute.residuals = FALSE,
    dims = 1:50,
    verbose = F)

predictions <- TransferData(
    anchorset = transfer_anchors, 
    refdata = referencehippocampus$subclass_label, 
    weight.reduction = combined.sct[['pca']],
    dims = 1:50,
    verbose = F)

combined.sct <- AddMetaData(
    object = combined.sct,
    metadata = predictions)
    
combined.sct$atlas_predictions = combined.sct$predicted.id

# Add cell cycle scores 
load("ref/mouse_cellcycle_genes.rda")
DefaultAssay(combined.sct) = "SCT"
combined.sct<- CellCycleScoring(combined.sct, s.features = m.s.genes, g2m.features = m.g2m.genes)

saveRDS(combined.sct,file="seurat/hippocampus_Parse_10x_integrated.rds")

sessionInfo()
