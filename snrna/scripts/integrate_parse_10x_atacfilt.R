library(Matrix)
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(viridis))
library(glmGamPoi)
library(RColorBrewer)
library(optparse)
library(stringr)
options(future.globals.maxSize = 18000 * 1024^2)


option_list = list(
  make_option(c("-t", "--tissue"), action="store", default=NA, type='character',
              help="hippocampus, cortex, heart, adrenal, or gastrocnemius"))

opt = parse_args(OptionParser(option_list=option_list))
tissue = opt$tissue


setwd("../../snrna/")


# also filter 10x for nuclei passing filters in snATAC
atac_metadata = read.csv("../snatac/ref/atac_metadata_all_tissues.csv")

tissue = tolower(tissue)
path = paste0("seurat/",str_to_lower(tissue),"_Parse_10x_integrated.rds")
combined.sct = readRDS(path)

obj_parse = subset(combined.sct, subset = technology == "Parse")
obj_parse_standard = subset(obj_parse, subset = depth1  == "shallow")
obj_parse_deep = subset(obj_parse, subset = depth1  == "deep")

obj_10x = subset(combined.sct, subset = technology == "10x")

obj_10x <- subset(obj_10x, subset = cellID %in% atac_metadata$cellID)


# SCT + CCA normalization and integration 
#Use pretty standard Seurat pipeline to perform SCT normalization and integration. Create list of the 3 Seurat objects, use additional package to make SCT go faster, and save pre-integrated data in `seurat` folder. Use Parse standard Seurat object as reference dataset because it contains all timepoints, while 10x data only contains 2 timepoints.

obj.list = list(obj_10x,obj_parse_standard,obj_parse_deep)

obj.list <- lapply(X = obj.list, FUN = SCTransform, method = "glmGamPoi", 
                   vars.to.regress = c("percent.mt","nFeature_RNA"), verbose = F)

saveRDS(obj.list,file=paste0("seurat/",str_to_lower(tissue),"_Parse_10x_to_integrate.rds"))

features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000, verbose = F)
obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = features, verbose = F)

names(obj.list) = c("10x","standard","deep")

reference_dataset <- which(names(obj.list) == "standard")

anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT", 
    anchor.features = features, reference = reference_dataset, verbose = F)
combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = F)

saveRDS(combined.sct,file=paste0("seurat/",str_to_lower(tissue),"_Parse_10x_integrated.rds"))

#Dimensionality reduction and clustering
#Standard Seurat processing with PCA, UMAP, SNN graph construction, and clustering. Use high clustering resolution to separate smaller subtypes.

DefaultAssay(combined.sct) = "integrated" # Make sure to cluster on the integrated assay

# PCA
combined.sct <- RunPCA(combined.sct, verbose = T, npcs = 50)

# UMAP and clustering
combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:30,verbose = F)
combined.sct <- FindNeighbors(combined.sct, reduction = "pca", dims = 1:30,verbose = F)
combined.sct <- FindClusters(combined.sct,resolution=1.6,verbose = F)

# Add cell cycle scores 
load("ref/mouse_cellcycle_genes.rda")
DefaultAssay(combined.sct) = "SCT"
combined.sct<- CellCycleScoring(combined.sct, s.features = m.s.genes, g2m.features = m.g2m.genes)

# SAVE
saveRDS(combined.sct,file=paste0("seurat/",str_to_lower(tissue),"_Parse_10x_integrated.rds"))
write.csv(combined.sct@meta.data,file=paste0("seurat/",str_to_lower(tissue),"_Parse_10x_integrated_metadata.csv"))

## Call cluster marker genes 
DefaultAssay(combined.sct)= "RNA"
Idents(combined.sct) = "seurat_clusters"
markers <- FindAllMarkers(combined.sct, 
                             only.pos = TRUE, 
                             min.pct = 0.25, 
                             logfc.threshold = 0.25, 
                             verbose = T)

write.table(markers,file=paste0("seurat/markers/",str_to_lower(tissue),"_cluster_marker_genes_only.pos_min.pct0.25_logfc.threshold0.25.tsv"),
            sep="\t",
            quote=F,
            row.names = F)

sessionInfo()
