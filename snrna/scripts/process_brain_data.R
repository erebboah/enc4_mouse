library(Matrix)
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(viridis))
library(glmGamPoi)
library(RColorBrewer)
options(future.globals.maxSize = 10000 * 1024^2)
setwd("/share/crsp/lab/seyedam/share/enc4_mouse/snrna/")

ref = readMM("ref/external_data/brain_atlas_1000_per_subtype_250734.mtx")

genes = read.delim("ref/external_data/brain_atlas_1000_per_subtype_250734_genes.tsv",header=F)

barcodes = read.delim("ref/external_data/brain_atlas_1000_per_subtype_250734_barcodes.tsv",header=F)

rownames(ref) = genes$V1
colnames(ref) = barcodes$V1

ref = ref[1:(nrow(ref)-129),]

meta = read.csv("ref/external_data/brain_atlas_metadata_1000_per_subtype_250734.csv")

ref_obj = CreateSeuratObject(counts = ref, min.cells = 0, min.features = 0)
ref_obj[["percent.mt"]] = PercentageFeatureSet(ref_obj, pattern = "^mt-")
ref_obj@meta.data = cbind(ref_obj@meta.data,meta)

ref_obj <- SCTransform(ref_obj, vars.to.regress = c("nFeature_RNA","percent.mt"),
                  method = "glmGamPoi", conserve.memory = TRUE)


ref_obj <- RunPCA(ref_obj, verbose = T, assay = "SCT")

saveRDS(ref_obj,file="ref/external_data/brain_atlas_subsampled_sct_pca.rds")
