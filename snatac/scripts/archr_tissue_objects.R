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

# adrenal gland

adr = proj$cellNames[proj_meta$tissue == "Adrenal"]

adr_proj = subsetArchRProject(
    ArchRProj = proj,
    outputDirectory = "ENC4_Mouse_Adrenal",
    cells = adr,
    dropCells = TRUE,
    force = TRUE)

## archr processing
adr_proj = addIterativeLSI(adr_proj, sampleCellsPre = 40000, varFeatures = 100000,
                       useMatrix = "TileMatrix",name = "IterativeLSI", force = TRUE)

adr_proj = addClusters(adr_proj, maxClusters = 40, 
                   reducedDims = "IterativeLSI",
                   force = TRUE)

adr_proj = addUMAP(adr_proj, reducedDims = "IterativeLSI", force = TRUE)

p1 <- plotEmbedding(ArchRProj = adr_proj, colorBy = "cellColData", name = "sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = adr_proj, colorBy = "cellColData", name = "tissue", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = adr_proj, colorBy = "cellColData", name = "timepoint", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = adr_proj, colorBy = "cellColData", name = "gen_celltype", embedding = "UMAP")
p5 <- plotEmbedding(ArchRProj = adr_proj, colorBy = "cellColData", name = "celltypes", embedding = "UMAP")
p6 <- plotEmbedding(ArchRProj = adr_proj, colorBy = "cellColData", name = "subtypes", embedding = "UMAP")

plotPDF(p1,p2,p3,p4,p5,p6, name = "Plot-UMAP-Sample-Clusters.pdf",
        ArchRProj = adr_proj, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(
  ArchRProj = adr_proj,
  outputDirectory = "ENC4_Mouse_Adrenal/"
)

print("Adrenal done")

# cortex

proj = loadArchRProject(path = "ENC4_Mouse/")
proj_meta = as.data.frame(proj@cellColData)

ctx = proj$cellNames[proj_meta$tissue == "Cortex"]

ctx_proj = subsetArchRProject(
    ArchRProj = proj,
    outputDirectory = "ENC4_Mouse_Cortex/",
    cells = ctx,
    dropCells = TRUE,
    force = TRUE)

## archr processing
ctx_proj = addIterativeLSI(ctx_proj, sampleCellsPre = 40000, varFeatures = 100000,
                       useMatrix = "TileMatrix",name = "IterativeLSI", force = TRUE)

ctx_proj = addClusters(ctx_proj, maxClusters = 40,
                   reducedDims = "IterativeLSI",
                   force = TRUE)

ctx_proj = addUMAP(ctx_proj, reducedDims = "IterativeLSI", force = TRUE)

p1 <- plotEmbedding(ArchRProj = ctx_proj, colorBy = "cellColData", name = "sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = ctx_proj, colorBy = "cellColData", name = "tissue", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = ctx_proj, colorBy = "cellColData", name = "timepoint", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = ctx_proj, colorBy = "cellColData", name = "gen_celltype", embedding = "UMAP")
p5 <- plotEmbedding(ArchRProj = ctx_proj, colorBy = "cellColData", name = "celltypes", embedding = "UMAP")
p6 <- plotEmbedding(ArchRProj = ctx_proj, colorBy = "cellColData", name = "subtypes", embedding = "UMAP")

plotPDF(p1,p2,p3,p4,p5,p6, name = "Plot-UMAP-Sample-Clusters.pdf",
        ArchRProj = ctx_proj, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(
  ArchRProj = ctx_proj,
  outputDirectory = "ENC4_Mouse_Cortex/"
)

print("Cortex done")

# hippocampus

proj = loadArchRProject(path = "ENC4_Mouse/")
proj_meta = as.data.frame(proj@cellColData)

hc = proj$cellNames[proj_meta$tissue == "Hippocampus"]

hc_proj = subsetArchRProject(
    ArchRProj = proj,
    outputDirectory = "ENC4_Mouse_Hippocampus/",
    cells = hc,
    dropCells = TRUE,
    force = TRUE)

## archr processing
hc_proj = addIterativeLSI(hc_proj, sampleCellsPre = 40000, varFeatures = 100000,
                       useMatrix = "TileMatrix",name = "IterativeLSI", force = TRUE)

hc_proj = addClusters(hc_proj, maxClusters = 40,
                   reducedDims = "IterativeLSI",
                   force = TRUE)

hc_proj = addUMAP(hc_proj, reducedDims = "IterativeLSI", force = TRUE)

p1 <- plotEmbedding(ArchRProj = hc_proj, colorBy = "cellColData", name = "sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = hc_proj, colorBy = "cellColData", name = "tissue", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = hc_proj, colorBy = "cellColData", name = "timepoint", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = hc_proj, colorBy = "cellColData", name = "gen_celltype", embedding = "UMAP")
p5 <- plotEmbedding(ArchRProj = hc_proj, colorBy = "cellColData", name = "celltypes", embedding = "UMAP")
p6 <- plotEmbedding(ArchRProj = hc_proj, colorBy = "cellColData", name = "subtypes", embedding = "UMAP")

plotPDF(p1,p2,p3,p4,p5,p6, name = "Plot-UMAP-Sample-Clusters.pdf",
        ArchRProj = hc_proj, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(
  ArchRProj = hc_proj,
  outputDirectory = "ENC4_Mouse_Hippocampus/"
)

print("Hippocampus done")

# heart

proj = loadArchRProject(path = "ENC4_Mouse/")
proj_meta = as.data.frame(proj@cellColData)

heart = proj$cellNames[proj_meta$tissue == "Heart"]

hrt_proj = subsetArchRProject(
    ArchRProj = proj,
    outputDirectory = "ENC4_Mouse_Heart/",
    cells = heart,
    dropCells = TRUE,
    force = TRUE)

## archr processing
hrt_proj = addIterativeLSI(hrt_proj, sampleCellsPre = 40000, varFeatures = 100000,
                       useMatrix = "TileMatrix",name = "IterativeLSI", force = TRUE)

hrt_proj = addClusters(hrt_proj, maxClusters = 40,
                   reducedDims = "IterativeLSI",
                   force = TRUE)

hrt_proj = addUMAP(hrt_proj, reducedDims = "IterativeLSI", force = TRUE)

p1 <- plotEmbedding(ArchRProj = hrt_proj, colorBy = "cellColData", name = "sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = hrt_proj, colorBy = "cellColData", name = "tissue", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = hrt_proj, colorBy = "cellColData", name = "timepoint", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = hrt_proj, colorBy = "cellColData", name = "gen_celltype", embedding = "UMAP")
p5 <- plotEmbedding(ArchRProj = hrt_proj, colorBy = "cellColData", name = "celltypes", embedding = "UMAP")
p6 <- plotEmbedding(ArchRProj = hrt_proj, colorBy = "cellColData", name = "subtypes", embedding = "UMAP")

plotPDF(p1,p2,p3,p4,p5,p6, name = "Plot-UMAP-Sample-Clusters.pdf",
        ArchRProj = hrt_proj, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(
  ArchRProj = hrt_proj,
  outputDirectory = "ENC4_Mouse_Heart/"
)

print("Heart done")

# gastrocnemius

proj = loadArchRProject(path = "ENC4_Mouse/")
proj_meta = as.data.frame(proj@cellColData)

gas = proj$cellNames[proj_meta$tissue == "Gastrocnemius"]

gas_proj = subsetArchRProject(
    ArchRProj = proj,
    outputDirectory = "ENC4_Mouse_Gastrocnemius/",
    cells = gas,
    dropCells = TRUE,
    force = TRUE)

## archr processing
gas_proj = addIterativeLSI(gas_proj, sampleCellsPre = 40000, varFeatures = 100000,
                       useMatrix = "TileMatrix",name = "IterativeLSI", force = TRUE)

gas_proj = addClusters(gas_proj, maxClusters = 40,
                   reducedDims = "IterativeLSI",
                   force = TRUE)

gas_proj = addUMAP(gas_proj, reducedDims = "IterativeLSI", force = TRUE)

p1 <- plotEmbedding(ArchRProj = gas_proj, colorBy = "cellColData", name = "sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = gas_proj, colorBy = "cellColData", name = "tissue", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = gas_proj, colorBy = "cellColData", name = "timepoint", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = gas_proj, colorBy = "cellColData", name = "gen_celltype", embedding = "UMAP")
p5 <- plotEmbedding(ArchRProj = gas_proj, colorBy = "cellColData", name = "celltypes", embedding = "UMAP")
p6 <- plotEmbedding(ArchRProj = gas_proj, colorBy = "cellColData", name = "subtypes", embedding = "UMAP")

plotPDF(p1,p2,p3,p4,p5,p6, name = "Plot-UMAP-Sample-Clusters.pdf",
        ArchRProj = gas_proj, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(
  ArchRProj = gas_proj,
  outputDirectory = "ENC4_Mouse_Gastrocnemius/"
)

print("Gastrocnemius done")
