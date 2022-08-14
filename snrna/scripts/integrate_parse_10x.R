library(Matrix)
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(viridis))
library(glmGamPoi)
library(RColorBrewer)
library(optparse)
library(stringr)
options(future.globals.maxSize = 10000 * 1024^2)


option_list = list(
  make_option(c("-t", "--tissue"), action="store", default=NA, type='character',
              help="hippocampus, cortex, heart, adrenal, or gastrocnemius"))

opt = parse_args(OptionParser(option_list=option_list))
tissue = opt$tissue


setwd("../../snrna/")
meta = read.delim("ref/enc4_mouse_snrna_metadata.tsv")

# Functions
# read in sparse matrix and assign row and column names
get_counts = function(batch){
    counts = readMM(paste0("scrublet/",batch,"_matrix.mtx"))
    barcodes = read.delim(paste0("scrublet/",batch,"_barcodes_scrublet.tsv"),header = F, 
                          col.names=c("cellID","doublet_scores","doublets"))
    
    features = read.delim(paste0("scrublet/",batch,"_genes.tsv"),header = F) 
    rownames(counts) = features$V1 
    colnames(counts) = barcodes$cellID
    out = counts
}

# read in associated metadata
get_metadata = function(batch){
    barcodes = read.delim(paste0("scrublet/",batch,"_barcodes_scrublet.tsv"),header = F, 
                          col.names=c("cellID","doublet_scores","doublets"))
    barcodes$library_accession = do.call("rbind", strsplit(barcodes$cellID, "[.]"))[,2]
    barcodes = left_join(barcodes,meta,by = "library_accession")
    out = barcodes
}

# merge the counts across experimental "batches"
# for example we sequenced 2 Parse "deep" libraries that should be combined into 1 counts matrix
# the technical batch effects between the standard and deep Parse libraries (and the Parse and 10x libraries) requires CCA integration

merge_counts = function(batches_list){
    matrix_list = list()
    for (i in 1:length(batches_list)){
        batch = batches_list[i]
        matrix_list[[i]] = get_counts(batch)
    }
    
    if (length(batches_list) < 2){
       matrix = matrix_list[[1]] 
       out = matrix
    } else {
        matrix = matrix_list[[1]]
        for (j in 2:length(batches_list)){
            matrix = RowMergeSparseMatrices(matrix,matrix_list[[j]])
        }
        out = matrix
    }
}

# merge the metadata across experimental "batches"
merge_metadata = function(batches_list){
    meta_list = list()
    for (i in 1:length(batches_list)){
        batch = batches_list[i]
        meta_list[[i]] = get_metadata(batch)
    }
    
    if (length(batches_list) < 2){
       meta = meta_list[[1]] 
       out = meta
    } else {
        meta = meta_list[[1]]
        for (j in 2:length(batches_list)){
            meta = rbind(meta,meta_list[[j]])
        }
        out = meta
    }
}

# make seurat object
seurat_obj = function(counts,metadata){
    obj = CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0)
    obj@meta.data = cbind(obj@meta.data,metadata)
    obj[["percent.mt"]] = PercentageFeatureSet(obj, pattern = "^mt-")
    obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa")
    out = obj
}

# Read in data
#Use functions defined above to create 3 Seurat objects: Parse standard, Parse deep, and 10x. Also make sure to get associated metadata, which includes the QC filter information.

meta = meta[meta$tissue == str_to_title(tissue),]

# get the experimental batches for 10x, Parse standard, and Parse deep
tenx_batches = unique(meta$experiment_batch[meta$technology == "10x"])

parse_standard_batches = unique(meta$experiment_batch[meta$technology == "Parse" & 
                                              meta$depth1 == "shallow"])

parse_deep_batches = unique(meta$experiment_batch[meta$technology == "Parse" & 
                                                  meta$depth1 == "deep"])


tenx_counts = merge_counts(tenx_batches)
tenx_meta = merge_metadata(tenx_batches)

parse_standard_counts = merge_counts(parse_standard_batches)
parse_standard_meta = merge_metadata(parse_standard_batches)

parse_deep_counts = merge_counts(parse_deep_batches)
parse_deep_meta = merge_metadata(parse_deep_batches)

# Make Seurat objects

obj_10x = seurat_obj(tenx_counts, 
                     tenx_meta)

obj_parse_standard = seurat_obj(parse_standard_counts, 
                                parse_standard_meta)

obj_parse_deep = seurat_obj(parse_deep_counts, 
                            parse_deep_meta)


# Filter
#Use QC information in metadata to filter by # UMIs and # genes detected per nucleus as well as doublet scores and percent mitochondrial gene expression. 

obj_10x <- subset(obj_10x, 
                  subset = nCount_RNA > unique(obj_10x$lower_nCount_RNA) & 
                  nCount_RNA < unique(obj_10x$upper_nCount_RNA)  & 
                  nFeature_RNA > unique(obj_10x$lower_nFeature_RNA) & 
                  doublet_scores < unique(obj_10x$upper_doublet_scores) & 
                  percent.mt < unique(obj_10x$upper_percent.mt)) 

obj_parse_standard <- subset(obj_parse_standard, 
                            subset = nCount_RNA > unique(obj_parse_standard$lower_nCount_RNA) & 
                            nCount_RNA < unique(obj_parse_standard$upper_nCount_RNA)  & 
                            nFeature_RNA > unique(obj_parse_standard$lower_nFeature_RNA) & 
                            doublet_scores < unique(obj_parse_standard$upper_doublet_scores) & 
                            percent.mt < unique(obj_parse_standard$upper_percent.mt))

obj_parse_deep <- subset(obj_parse_deep, 
                         subset = nCount_RNA > unique(obj_parse_deep$lower_nCount_RNA) & 
                         nCount_RNA < unique(obj_parse_deep$upper_nCount_RNA)  & 
                         nFeature_RNA > unique(obj_parse_deep$lower_nFeature_RNA) & 
                         doublet_scores < unique(obj_parse_deep$upper_doublet_scores) & 
                         percent.mt < unique(obj_parse_deep$upper_percent.mt))

             

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

## Call cluster marker genes 
DefaultAssay(combined.sct)= "SCT"
Idents(combined.sct) = "seurat_clusters"
markers <- FindAllMarkers(combined.sct, 
                             test.use = "MAST",
                             only.pos = TRUE, 
                             min.pct = 0.25, 
                             logfc.threshold = 0.25, 
                             verbose = T)

write.table(markers,file=paste0("seurat/",str_to_lower(tissue),"_cluster_marker_genes_mast_only.pos_min.pct0.25_logfc.threshold0.25.tsv",
            sep="\t",
            quote=F)


# SAVE
saveRDS(combined.sct,file=paste0("seurat/",str_to_lower(tissue),"_Parse_10x_integrated.rds"))

sessionInfo()
