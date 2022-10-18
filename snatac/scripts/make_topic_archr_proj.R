suppressPackageStartupMessages(library(ArchR))
library(parallel)
library(viridis)
suppressPackageStartupMessages(library(Seurat))
library(reshape2)
library(tidyverse)
addArchRThreads(threads = 8) 
addArchRGenome("mm10")
library(optparse)
set.seed(1234)

option_list = list(
  make_option(c("-t", "--tissue"), action="store", default=NA, type='character',
              help="hippocampus, cortex, heart, adrenal, or gastrocnemius"))

opt = parse_args(OptionParser(option_list=option_list))
tissue = opt$tissue

setwd("../archr/")

path = paste0("ENC4_Mouse_",str_to_title(tissue),"_Topics")
setwd(path)

# read in the topic fragment files made in step 2
inputFiles <- list.files(pattern = "\\.tsv.gz$")

sampleNames = do.call("rbind", strsplit(as.character(inputFiles), "_fragments"))[,1]

# create arrow files
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
  outputDirectory = path
)

saveArchRProject(ArchRProj = topic_proj, outputDirectory = path)

# add metadata
tissue_path = paste0("ENC4_Mouse_",str_to_title(tissue))
proj = loadArchRProject(path = paste0("../",tissue_path))

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
                                "Clusters",
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
topic_proj$orig_archr_clusters = meta$Clusters

# call peaks
topic_proj <- addGroupCoverages(ArchRProj = topic_proj, groupBy = "Sample")
pathToMacs2 <- findMacs2()
topic_proj <- addReproduciblePeakSet(
    ArchRProj = topic_proj, 
    groupBy = "Sample", 
    pathToMacs2 = pathToMacs2
)

# add peak matrix
topic_proj <- addPeakMatrix(topic_proj)

# save
saveArchRProject(ArchRProj = topic_proj, outputDirectory = path)

# calculate differential peak accessibility between Samples (topic clusters)
markersPeaks <- getMarkerFeatures(
    ArchRProj = topic_proj, 
    useMatrix = "PeakMatrix", 
    groupBy = "Sample",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# motif enrichment in differentially accessible peaks
topic_proj <- addMotifAnnotations(ArchRProj = topic_proj, 
                                  motifSet = "vierstra",
				  collection="archetype", 
                                  name = "Motif")

enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = topic_proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5" # can try adjusting this
  )

# plot preliminary results
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 20, transpose = TRUE)
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", 
        width = 8, height =16, ArchRProj = topic_proj, addDOC = FALSE)
