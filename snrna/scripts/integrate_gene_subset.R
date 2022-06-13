library(Matrix)
library(stringr)
library(Seurat)

setwd("/share/crsp/lab/seyedam/share/enc4_mouse/snrna/")
tfs = read.delim("ref/TF_mouse_GRCm39.txt")
tf_include = unique(tfs$Gene.name)
length(tf_include)

chrom = read.delim("ref/histone_regulators.tsv")
dna_include = unique(chrom$gene_name)
length(dna_include)

mirna = read.delim("ref/mirna_host_correlations.tsv")
hgs_include = unique(c(mirna$hg[mirna$spearman_correlation>0.3],
                   mirna$hg[grep("^M.+hg$",mirna$hg)]))
length(hgs_include)

genes_include = unique(c(tf_include,dna_include,hgs_include))
length(genes_include)

tissues = c("Adrenal","Gastrocnemius","Heart","Cortex","Hippocampus")

tissue = tissues[1]
obj = readRDS(paste0("seurat/",tissue,"_seurat.rds"))
counts = obj@assays$RNA@counts
counts = counts[rownames(counts) %in% genes_include,]
meta = obj@meta.data

for (j in 2:length(tissues)){
    tissue = tissues[j]
    obj = readRDS(paste0("seurat/",tissue,"_Parse_10x_integrated.rds"))
    
    raw_counts = obj@assays$RNA@counts # get raw counts
    raw_counts = raw_counts[rownames(raw_counts) %in% genes_include,]
    metadata = obj@meta.data

    # combine all tissues
    counts = RowMergeSparseMatrices(counts,raw_counts)
    meta = rbind(meta,metadata) 
} 


all = CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0)
all@meta.data = cbind(all@meta.data,meta[,!(colnames(meta) %in% c("X","nCount_RNA","nFeature_RNA","orig.ident"))])# dupe columns get generated


all_list = SplitObject(all, split.by = "tissue") 
all_list <- lapply(X = all_list, FUN = SCTransform, vars.to.regress = c("nFeature_RNA"))
features <- SelectIntegrationFeatures(object.list = all_list, nfeatures = 3000)
all_list <- PrepSCTIntegration(object.list = all_list, anchor.features = features)
all_list <- lapply(X = all_list, FUN = RunPCA, features = features)
anchors <- FindIntegrationAnchors(object.list = all_list, normalization.method = "SCT",
    anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)

combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30)
combined.sct <- RunPCA(combined.sct, verbose = FALSE)
combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:30)
combined.sct <- FindNeighbors(combined.sct, reduction = "pca", dims = 1:30)
combined.sct <- FindClusters(combined.sct)            
saveRDS(combined.sct,file=paste0(paste0("seurat/all_tissues_tfs_mirhgs_chromreg_Parse_10x_integrated.rds")))    

