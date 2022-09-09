suppressPackageStartupMessages(library(ArchR))
library(parallel)
library(viridis)
suppressPackageStartupMessages(library(Seurat))
set.seed(1234)
library(reshape2)
library(dplyr)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")

setwd("/share/crsp/lab/seyedam/share/enc4_mouse/snatac/topics/")

# read in the topic fragment files made in step 2
inputFiles <- list.files(pattern = "\\.tsv.gz$")

sampleNames = do.call("rbind", strsplit(as.character(inputFiles), "_fragments"))[,1]
sampleNames = do.call("rbind", strsplit(as.character(sampleNames), "topic_"))[,2]

# creat arrow files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = sampleNames,
  offsetPlus = 0,
  offsetMinus = 0,
  minTSS = 4,
  minFrags = 1000, 
  excludeChr = c("chrM"),
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

# make project
topic_proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "ENC4_Mouse_Topics_Heart"
)

# add metadata
setwd("/share/crsp/lab/seyedam/share/enc4_mouse/snatac/archr/")

proj = loadArchRProject(path = "ENC4_Mouse_heart/")

celltype_meta = getCellColData(proj)
celltype_meta$cellNames = rownames(celltype_meta)
celltype_meta = celltype_meta[,
                              c("technology",
                                "species",
                                "timepoint",
                                "sex",
                                "rep",
                                "tissue",
                               "sample",
                               "gen_celltype",
                               "celltypes",
                               "subtypes",
                               "rna_bc",
                               "cellID",
                               "cellNames")]
celltype_meta = as.data.frame(celltype_meta)

meta = getCellColData(topic_proj)
meta$origcellNames =rownames(meta)
meta$cellNames =rownames(meta)
meta$cellNames = do.call("rbind", strsplit(as.character(meta$cellNames), "#"))[,2]
meta$cellNames = do.call("rbind", strsplit(as.character(meta$cellNames), ":"))[,2]
meta$cellNames = gsub("[.]","#",meta$cellNames)
meta = as.data.frame(meta)

meta = left_join(meta,celltype_meta)

topic_proj$rna_bc = meta$rna_bc
topic_proj$cellID = meta$cellID
topic_proj$technology= meta$technology
topic_proj$species= meta$species
topic_proj$timepoint= meta$timepoint
topic_proj$sex= meta$sex
topic_proj$rep= meta$rep
topic_proj$tissue = meta$tissue
topic_proj$sample = meta$sample
topic_proj$gen_celltype = meta$gen_celltype
topic_proj$celltypes = meta$celltypes
topic_proj$subtypes = meta$subtypes

setwd("/share/crsp/lab/seyedam/share/enc4_mouse/snatac/topics/")

saveArchRProject(ArchRProj = topic_proj, outputDirectory = "ENC4_Mouse_Topics_Heart")

