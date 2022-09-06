library(Matrix)
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(viridis))
library(glmGamPoi)
library(RColorBrewer)
options(future.globals.maxSize = 10000 * 1024^2)
setwd("../../snrna/")

# Predict cell types from references

# Stressed heart (https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8810/files/)

counts = read.delim("ref/external_data/full_count_matrix.tsv")

ref = CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0)
ref[["percent.mt"]] = PercentageFeatureSet(ref, pattern = "^mt-")
ref <- SCTransform(ref, vars.to.regress = c("nFeature_RNA","percent.mt"))
ref <- RunPCA(ref, verbose = T, assay = "SCT")
ref <- RunUMAP(ref, reduction = "pca", dims = 1:30,verbose = F)
ref <- FindNeighbors(ref, reduction = "pca", dims = 1:30,verbose = F)
ref <- FindClusters(ref,resolution=0.25,verbose = F)
genes = c("Ms4a1", # b cells
          "Ttn",  # cardiomyocytes
          "Cd209a",# dc-like
          "Npr3","Pecam1", # endocardial; endothelial just pecam
          "Rgcc", # endothelial
          "Msln", # epicardial
          "Pdgfra", # fibroblasts
          "Wif1", # fibroblasts wif1+
          "S100a9", # granulocytes
          "Lyve1", # lymphatic ec
          "Fcgr1", # macrophages
          "Plp1", # schwann cells
          "Acta2","Pdgfrb", # smooth muscle
          "Ncr1","Cd3e" # t cells
         )

# dot plot of marker genes
DefaultAssay(ref) = "SCT"

# https://www.ahajournals.org/doi/10.1161/CIRCULATIONAHA.119.045115
pdf(file="plots/heart/annotation/ref_marker_cluster_dotplot.pdf",
    width = 10, height = 5)
DotPlot(ref, features = genes)+ 
theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
dev.off()

ref$celltypes = ref$seurat_clusters
ref$celltypes = gsub("\\<0\\>","Fibroblasts",ref$celltypes)
ref$celltypes = gsub("\\<1\\>","Macrophages",ref$celltypes)
ref$celltypes = gsub("\\<2\\>","Endothelial",ref$celltypes)
ref$celltypes = gsub("\\<3\\>","Fibroblasts_Wif1+",ref$celltypes)
ref$celltypes = gsub("\\<4\\>","Fibroblasts",ref$celltypes)
ref$celltypes = gsub("\\<5\\>","Smooth_muscle",ref$celltypes)
ref$celltypes = gsub("\\<6\\>","Endocardial",ref$celltypes)
ref$celltypes = gsub("\\<7\\>","Fibroblasts",ref$celltypes)
ref$celltypes = gsub("\\<8\\>","DC-like_cells",ref$celltypes)
ref$celltypes = gsub("\\<9\\>","Cardiomyocytes",ref$celltypes)
ref$celltypes = gsub("\\<10\\>","Granulocytes",ref$celltypes)
ref$celltypes = gsub("\\<11\\>","B_cells",ref$celltypes)
ref$celltypes = gsub("\\<12\\>","T_cells",ref$celltypes)
ref$celltypes = gsub("\\<13\\>","Lymphatic_ECs",ref$celltypes)
ref$celltypes = gsub("\\<14\\>","Schwann_cells",ref$celltypes)
ref$celltypes = gsub("\\<15\\>","Epicardial",ref$celltypes)
ref$celltypes = gsub("\\<16\\>","Fibroblasts",ref$celltypes)

pdf(file="plots/heart/annotation/ref_marker_celltype_dotplot.pdf",
    width = 10, height = 5)
DefaultAssay(ref) = "SCT"
Idents(ref) ="celltypes"
Idents(ref) = factor(Idents(ref), levels = sort(as.character(unique(Idents(ref)))))

DotPlot(ref, features = genes)+ 
theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
dev.off()

saveRDS(ref,file="ref/external_data/E-MTAB-8810_processed.rds")

combined.sct = readRDS("seurat/heart_Parse_10x_integrated.rds")

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

# Add cell cycle scores 
load("ref/mouse_cellcycle_genes.rda")
DefaultAssay(combined.sct) = "SCT"
combined.sct<- CellCycleScoring(combined.sct, s.features = m.s.genes, g2m.features = m.g2m.genes)

saveRDS(combined.sct,file="seurat/heart_Parse_10x_integrated.rds")
rm(ref) 

# Human atlas converted to mouse

library(SeuratDisk)
Convert("ref/external_data/global_raw.h5ad", dest = "h5seurat", overwrite = TRUE)
ref_human <- LoadH5Seurat("ref/external_data/global_raw.h5seurat")
ref_counts = ref_human@assays$RNA@counts

# Keep getting biomaRt server error, but this would be the correct way--
# Basic function to convert human to mouse gene names
#convertHumanGeneList <- function(x){
#require("biomaRt")
#human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
#genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
#humanx <- genes$V2
#return(humanx)
#}
#genes <- convertHumanGeneList(rownames(ref_counts))
#ref_counts = ref_counts[rownames(ref_counts) %in% genes$HGNC.symbol,]
#genes=genes[match(rownames(ref_counts),genes$HGNC.symbol),]
#rownames(ref_counts) = genes$MGI.symbol

# simply convert to mouse gene format
rownames(ref_counts) = str_to_title(rownames(ref_counts))

ref_meta = ref_human@meta.data
rm(ref_human)  
                         
ref_atlas = CreateSeuratObject(counts = ref_counts, min.cells = 0, min.features = 0)
ref_atlas@meta.data = cbind(ref_atlas@meta.data,ref_meta)
ref_atlas[["percent.mt"]] = PercentageFeatureSet(ref_atlas, pattern = "^mt-")
ref_atlas <- SCTransform(ref_atlas, method = "glmGamPoi", conserve.memory = TRUE, vars.to.regress = c("nFeature_RNA","percent.mt"))
ref_atlas <- RunPCA(ref_atlas, verbose = T, assay = "SCT")

DefaultAssay(combined.sct) <- "SCT"
DefaultAssay(ref_atlas) <- "SCT"

saveRDS(ref_atlas,file="ref/external_data/heart_atlas_processed.rds")

transfer_anchors <- FindTransferAnchors(
    reference = ref_atlas,
    query = combined.sct,
    reference.assay = "SCT",
    normalization.method = "SCT",
    reference.reduction = "pca",
    recompute.residuals = FALSE,
    dims = 1:50,
    verbose = F)

predictions <- TransferData(
    anchorset = transfer_anchors, 
    refdata = ref_atlas$cell_type, 
    weight.reduction = combined.sct[['pca']],
    dims = 1:50,
    verbose = F)

combined.sct <- AddMetaData(
    object = combined.sct,
    metadata = predictions)
    
combined.sct$atlas_predictions = combined.sct$predicted.id                        

saveRDS(combined.sct,file="seurat/heart_Parse_10x_integrated.rds")

sessionInfo()
