library(Seurat)
suppressPackageStartupMessages(library(ArchR))
set.seed(1234)
library(parallel)
addArchRThreads(16)
addArchRGenome("mm10")
library(tidyverse)

setwd("../archr/")

proj = loadArchRProject(path = "ENC4_Mouse/")
proj_meta = as.data.frame(proj@cellColData)

metadata = read.delim("../ref/enc4_mouse_snatac_metadata.tsv")
all_tissues_10x_rna_metadata = read.csv("../../snrna/seurat/all_tissues_Parse_10x_TFs_mirhgs_chromreg_metadata.csv",row.names = "X")

print("how many ATAC cells passed RNA filters?")
table(proj_meta$cellID %in% all_tissues_10x_rna_metadata$cellID)

# filter cells from ATAC that are not in RNA
# also QC filtere
atac_cells_filt = proj$cellNames[proj_meta$cellID %in% all_tissues_10x_rna_metadata$cellID]

proj = subsetArchRProject(
    ArchRProj = proj,
    outputDirectory = "ENC4_Mouse/",
    cells = atac_cells_filt,
    dropCells = TRUE,
    force = TRUE)

# archr processing
proj = addIterativeLSI(proj, sampleCellsPre = 40000, varFeatures = 100000,
                       useMatrix = "TileMatrix",name = "IterativeLSI", force = T)

proj = addClusters(proj, maxClusters = 40, 
                   reducedDims = "IterativeLSI",
                   force = TRUE)

proj = addUMAP(proj, reducedDims = "IterativeLSI", force = T)

# add metadata from RNA cell type annotation
proj_meta = as.data.frame(proj@cellColData)
proj_meta = as.data.frame(proj_meta[,c("cellID")])
colnames(proj_meta) = "cellID"
proj_meta = dplyr::left_join(proj_meta, all_tissues_10x_rna_metadata, "cellID") # merge by cellID
table(proj_meta$cellID == proj$cellID) # sanity check

# add RNA metadata to ATAC
proj$technology= proj_meta$technology
proj$species= proj_meta$species
proj$timepoint= proj_meta$timepoint
proj$sex= proj_meta$sex
proj$rep= proj_meta$rep
proj$tissue = proj_meta$tissue
proj$sample = proj_meta$sample
proj$gen_celltype = proj_meta$gen_celltype
proj$celltypes = proj_meta$celltypes
proj$subtypes = proj_meta$subtypes

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "tissue", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "timepoint", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "gen_celltype", embedding = "UMAP")
p5 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "celltypes", embedding = "UMAP")
p6 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "subtypes", embedding = "UMAP")

plotPDF(p1,p2,p3,p4,p5,p6, name = "Plot-UMAP-Sample-Clusters.pdf", 
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(
  ArchRProj = proj,
  outputDirectory = "ENC4_Mouse/"
)

