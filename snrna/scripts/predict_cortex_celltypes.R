library(Matriix)
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(viridis))
library(glmGamPoi)
library(RColorBrewer)
options(future.globals.maxSize = 10000 * 1024^2)
setwd("../../snrna/")

# Predict cell types from references

ref = readRDS("ref/external_data/brain_atlas_subsampled_sct_pca.rds")
combined.sct = readRDS("seurat/cortex_Parse_10x_integrated.rds")

DefaultAssay(combined.sct) <- "SCT"
DefaultAssay(ref) <- "SCT"

transfer_anchors <- FindTransferAnchors(
    reference = ref,
    query = combined.sct,
    reference.assay = "SCT",
    normalization.method = "SCT",
    reference.reduction = "pca",
    recompute.residuals = FALSE,
    dims = 1:50,
    verbose = F)

predictions <- TransferData(
    anchorset = transfer_anchors, 
    refdata = ref$subclass_label, 
    weight.reduction = combined.sct[['pca']],
    dims = 1:50,
    verbose = F)

combined.sct <- AddMetaData(
    object = combined.sct,
    metadata = predictions)
    
combined.sct$atlas_predictions = combined.sct$predicted.id

saveRDS(combined.sct,file="seurat/cortex_Parse_10x_integrated.rds")

sessionInfo()
