library(Seurat)
suppressPackageStartupMessages(library(ArchR))
set.seed(1234)
library(parallel)
addArchRThreads(8)
addArchRGenome("mm10")

setwd("/share/crsp/lab/seyedam/share/enc4_mouse/snatac/archr/")
metadata = read.delim("../ref/enc4_mouse_snatac_metadata.tsv")

inputFiles = paste0("../fragments/",metadata$file_accession,
                       "/encode_scatac_dcc_2/results/",
                        metadata$experiment_accession,"-1/fragments/fragments.tsv.gz")

names(inputFiles) = metadata$sample

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
  #force = T
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "ENC4_Mouse"
)

proj$sample = do.call("rbind", strsplit(as.character(proj$cellNames), "#"))[,1]
proj_meta = as.data.frame(proj@cellColData)
proj_meta = plyr::join(proj_meta, metadata)
proj$library_accession= proj_meta$library_accession
proj$rna_library_accession= proj_meta$rna_library_accession
proj$rna_experiment_accession= proj_meta$rna_experiment_accession
proj$experiment_accession= proj_meta$experiment_accession
proj$file_accession= proj_meta$file_accession
proj$multiome_experiment_accession= proj_meta$multiome_experiment_accession
proj$timepoint= proj_meta$timepoint
proj$sex= proj_meta$sex
proj$rep= proj_meta$rep
proj$tissue = proj_meta$tissue

# read in 10x barcode sequences
atac_bcs = read.delim("../ref/atac_737K-arc-v1.txt",header = F)
rna_bcs = read.delim("../ref/gene_exp_737K-arc-v1.txt",header = F)
bcs = cbind(atac_bcs,rna_bcs)
colnames(bcs) = c("atac_bc","rna_bc")

proj$cellbarcode = do.call("rbind", strsplit(as.character(proj$cellNames), "#"))[,2]
proj$atac_bc = toupper(spgs::reverseComplement(proj$cellbarcode))

bcs = bcs[bcs$atac_bc %in% proj$atac_bc,]
bcs = bcs[match(proj$atac_bc, bcs$atac_bc),]

proj$rna_bc = paste0(bcs$rna_bc,".",proj$rna_library_accession)

# filter archr nuclei
all = readRDS("/share/crsp/lab/seyedam/share/Heidi_Liz/all_mouse/all_tissues_Parse_10x_integrated_TFs_mirhgs.rds")

cells = all$cellID
meta = all@meta.data[,c("cellID","celltypes","gen_celltype","subtypes")]
rna_meta = meta

atac_cells_filt = proj$cellNames[proj$rna_bc %in% rna_meta$cellID]

proj = subsetArchRProject(
    ArchRProj = proj,
    cells = atac_cells_filt,
    dropCells = TRUE, #  drop cells that are not in ArchRProject from corresponding Arrow Files 
    force = TRUE)

# add cell type metadata
rna_meta = rna_meta[rna_meta$cellID %in% proj$rna_bc,]
rna_meta = rna_meta[match(proj$rna_bc,rna_meta$cellID),]

table(proj$rna_bc == rna_meta$cellID)
proj$subtypes = rna_meta$subtypes
proj$celltypes = rna_meta$celltypes
proj$gen_celltype = rna_meta$gen_celltype

# add doublet scores and filter
proj = addDoubletScores(proj, k = 10, knnMethod = "LSI")
proj <- filterDoublets(proj)

# arhr processing
proj = addIterativeLSI(proj, sampleCellsPre = 40000, varFeatures = 100000,
                       useMatrix = "TileMatrix",name = "IterativeLSI", force = T)

proj = addClusters(proj, maxClusters = 40, 
                   reducedDims = "IterativeLSI",
                   force = TRUE)

proj = addUMAP(proj, reducedDims = "IterativeLSI", force = T)

# plots
p1 <- plotEmbedding(
  proj, 
  colorBy = "cellColData", 
  name = "sample"
)

p2 <- plotEmbedding(
  proj, 
  colorBy = "cellColData", 
  name = "sex"
)

p3<- plotEmbedding(
  proj, 
  colorBy = "cellColData", 
  name = "timepoint"
)

p4 <- plotEmbedding(
  proj, 
  colorBy = "cellColData", 
  name = "Clusters"
)

p5 <- plotEmbedding(
  proj, 
  colorBy = "cellColData", 
  name = "subtypes"
)

p6 <- plotEmbedding(
  proj, 
  colorBy = "cellColData", 
  name = "celltypes"
)

p7 <- plotEmbedding(
  proj, 
  colorBy = "cellColData", 
  name = "gen_celltype"
)

pdf(file="ENC4_Mouse/atac_all_tissue_sample_cluster_umaps.pdf",width=10,height=5)
p1
p2
p3
p4
dev.off()

pdf(file="ENC4_Mouse/atac_all_tissue_celltype_umaps.pdf",width=10,height=10)
p7
p6
p5
dev.off()

options(repr.plot.width=18,repr.plot.height=8)

p1 <- plotEmbedding(
  proj, 
  colorBy = "cellColData", 
  name = "TSSEnrichment"
)

p2 <- plotEmbedding(
  proj, 
  colorBy = "cellColData", 
  name = "nFrags"
)

p3 <- plotEmbedding(
  proj, 
  colorBy = "cellColData", 
  name = "DoubletScore"
)

p4 <- plotEmbedding(
  proj, 
  colorBy = "cellColData", 
  name = "DoubletEnrichment"
)

pdf(file="ENC4_Mouse/atac_all_tissue_qc_umaps.pdf",width=10,height=5)
p1
p2
p3
p4
dev.off()

setwd("/share/crsp/lab/seyedam/share/enc4_mouse/snatac/archr/")

saveArchRProject(
  ArchRProj = proj,
  outputDirectory = "ENC4_Mouse/"
)
