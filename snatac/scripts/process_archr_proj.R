library(Seurat)
suppressPackageStartupMessages(library(ArchR))
set.seed(1234)
library(parallel)
addArchRThreads(16)
addArchRGenome("mm10")
library(tidyverse)

setwd("/share/crsp/lab/seyedam/share/enc4_mouse/snatac/archr/")

proj = loadArchRProject(path = "ENC4_Mouse/")

metadata = read.delim("../ref/enc4_mouse_snatac_metadata.tsv")
all_tissues_10x_rna_metadata = read.delim("/share/crsp/lab/seyedam/share/Heidi_Liz/all_mouse/all_tissues_10x_TFs_mirhgs_chromreg_metadata.tsv")

proj_meta = as.data.frame(proj@cellColData)

print("how many ATAC cells passed RNA filters?")
table(proj_meta$cellID %in% all_tissues_10x_rna_metadata$cellID)

# filter cells from ATAC that are not in RNA
# also QC filter
atac_cells_filt = proj$cellNames[proj_meta$cellID %in% all_tissues_10x_rna_metadata$cellID & proj_meta$TSSEnrichment > 4 & proj_meta$nFrags > 1000]

proj = subsetArchRProject(
    ArchRProj = proj,
    cells = atac_cells_filt,
    dropCells = TRUE, #  drop cells that are not in ArchRProject from corresponding Arrow Files 
    force = TRUE)

# archr processing
proj = addIterativeLSI(proj, sampleCellsPre = 40000, varFeatures = 100000,
                       useMatrix = "TileMatrix",name = "IterativeLSI", force = T)

proj = addClusters(proj, maxClusters = 40, 
                   reducedDims = "IterativeLSI",
                   force = TRUE)

proj = addUMAP(proj, reducedDims = "IterativeLSI", force = T)

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "tissue", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "timepoint", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "gen_celltype", embedding = "UMAP")
p5 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "celltypes", embedding = "UMAP")
p6 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "subtypes", embedding = "UMAP")

plotPDF(p1,p2,p3,p4,p5,p6 name = "Plot-UMAP-Sample-Clusters.pdf", 
        ArchRProj = topic_proj, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(
  ArchRProj = proj,
  outputDirectory = "ENC4_Mouse/"
)

