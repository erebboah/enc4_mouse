library(Seurat)
suppressPackageStartupMessages(library(ArchR))
set.seed(1234)
library(parallel)
addArchRThreads(8)
addArchRGenome("mm10")

setwd("/share/crsp/lab/seyedam/share/enc4_mouse/snatac/archr/")
metadata = read.delim("../ref/enc4_mouse_snatac_metadata.tsv")
all_tissues_10x_rna_metadata = read.delim("/share/crsp/lab/seyedam/share/Heidi_Liz/all_mouse/all_tissues_10x_TFs_mirhgs_chromreg_metadata.tsv")

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
  addGeneScoreMat = TRUE
)

# make archr project
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "ENC4_Mouse")

saveArchRProject(
  ArchRProj = proj,
  outputDirectory = "ENC4_Mouse/"
)

# add doublet scores and filter
proj = addDoubletScores(proj, k = 10, knnMethod = "LSI")
proj <- filterDoublets(proj)

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
proj_meta = dplyr::left_join(proj_meta, all_tissues_10x_rna_metadata) # merge by cellID and library accession

print("checking that new metadata is ordered correctly:")
table(proj_meta$cellNames == proj$cellNames)
print("how many ATAC cells passed RNA filters?")
table(proj_meta$cellID %in% all_tissues_10x_rna_metadata$cellID)

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

# filter cells from ATAC that are not in RNA
atac_cells_filt = proj$cellNames[proj_meta$cellID %in% all_tissues_10x_rna_metadata$cellID]

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

saveArchRProject(
  ArchRProj = proj,
  outputDirectory = "ENC4_Mouse/"
)

