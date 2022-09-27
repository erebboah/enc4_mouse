library(Matrix)
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(viridis))
library(glmGamPoi)
library(RColorBrewer)
options(future.globals.maxSize = 10000 * 1024^2)
setwd("../../snrna/")

# Predict cell types from references

combined.sct = readRDS("seurat/gastrocnemius_Parse_10x_integrated.rds")

# Adjusted UMAP function to increase cluster separation
DefaultAssay(combined.sct) <- "integrated"
combined.sct <- RunUMAP(combined.sct, reduction = "pca", spread = 2, dims = 1:30,verbose = F)
combined.sct <- FindNeighbors(combined.sct, reduction = "pca", dims = 1:30,verbose = F)
combined.sct <- FindClusters(combined.sct,resolution=2,verbose = F)


p10 = readRDS("ref/external_data/p10.rds")
p21 = readRDS("ref/external_data/p21.rds")
mo5 = readRDS("ref/external_data/5mo.rds")

ref = merge(p10, y = c(p21,mo5), add.cell.ids = c("PND10", "PND21", "PNM05"), project = "TA")
ref$celltypes = Idents(ref)
ref[["percent.mt"]] = PercentageFeatureSet(ref, pattern = "^mt-")
ref <- SCTransform(ref,  method = "glmGamPoi", 
                   vars.to.regress = c("percent.mt","nFeature_RNA"), verbose = F)
ref<- RunPCA(ref, verbose = T, npcs = 50)


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
    refdata = ref$celltypes, 
    weight.reduction = combined.sct[['pca']],
    dims = 1:50,
    verbose = F)

combined.sct <- AddMetaData(
    object = combined.sct,
    metadata = predictions)
    
combined.sct$predictions = combined.sct$predicted.id

saveRDS(combined.sct,file="seurat/gastrocnemius_Parse_10x_integrated.rds")

sessionInfo()
