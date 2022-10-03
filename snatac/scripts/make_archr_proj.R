library(Seurat)
suppressPackageStartupMessages(library(ArchR))
set.seed(1234)
library(parallel)
addArchRThreads(16)
addArchRGenome("mm10")
library(tidyverse)

setwd("../archr/")
metadata = read.delim("../ref/enc4_mouse_snatac_metadata.tsv")

# read in fragment files from ENCODE portal
inputFiles = paste0("../fragments/",metadata$file_accession,
                       "/encode_scatac_dcc_2/results/",
                        metadata$experiment_accession,"-1/fragments/fragments.tsv.gz")

# set names to be the RNA library ID to match the RNA barcodes later
names(inputFiles) = metadata$rna_library_accession

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  offsetPlus = 0,
  offsetMinus = 0,
  minTSS = 4,
  minFrags = 1000, 
  excludeChr = c("chrM"),
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force = TRUE
)

# make archr project
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "ENC4_Mouse")

saveArchRProject(
  ArchRProj = proj,
  outputDirectory = "ENC4_Mouse/"
)

# match ATAC barcodes to RNA barcodes
proj_meta = as.data.frame(proj@cellColData)
proj_meta$cellNames = rownames(proj_meta)
proj_meta$library_accession = proj_meta$Sample

# read in 10x barcode sequences
atac_bcs = read.delim("../ref/atac_737K-arc-v1.txt",header = F)
rna_bcs = read.delim("../ref/gene_exp_737K-arc-v1.txt",header = F)
bcs = cbind(atac_bcs,rna_bcs)
colnames(bcs) = c("atac_bc","rna_bc")

# match RNA and ATAC cell barcodes
proj_meta$atac_bc = do.call("rbind", strsplit(as.character(proj_meta$cellNames), "#"))[,2]
proj_meta$atac_bc = toupper(spgs::reverseComplement(proj_meta$atac_bc))

proj_meta = dplyr::left_join(proj_meta, bcs) # merge by atac_bc
proj_meta$cellID = paste0(proj_meta$rna_bc,".",proj_meta$library_accession) # re-create RNA cell IDs

print("checking that new metadata is ordered correctly:")
table(proj_meta$cellNames == proj$cellNames)
proj$original_cellNames = proj_meta$cellNames
proj$original_atac_bc= do.call("rbind", strsplit(proj$original_cellNames, "#"))[,2]
proj$atac_bc = proj_meta$atac_bc
proj$cellID = proj_meta$cellID

write.csv(as.data.frame(proj@cellColData),file="../ref/atac_unfiltered_metadata_all_tissues.csv")

# add doublet scores and filter
proj = addDoubletScores(proj, k = 10, knnMethod = "LSI")
proj <- filterDoublets(proj)

write.csv(as.data.frame(proj@cellColData),file="../ref/atac_metadata_all_tissues.csv")

saveArchRProject(
  ArchRProj = proj,
  outputDirectory = "ENC4_Mouse/"
)
