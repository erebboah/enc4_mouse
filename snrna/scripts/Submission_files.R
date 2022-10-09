# Create and save files for Synapse submissions
# - EAID_000084: Adrenal gland RNA + Multiome
# - EAID_000085: Left cortex RNA + Multiome
# - EAID_000086: Hippocampus RNA + Multiome
# - EAID_000087: Heart RNA + Multiome
# - EAID_000088: Gastrocnemius RNA + Multiome

# Guidelines: https://docs.google.com/document/d/1NAuihlO3Qd4y588Fpv1iWTyO2jxo2u8M9U-vfwnGbmE/edit

library(tidyverse)
library(Matrix)
library(Seurat)
library(ArchR)
library(data.table)

setwd("/share/crsp/lab/seyedam/share/enc4_mouse/")

system("mkdir synapse")
system("mkdir synapse/header")
system("mkdir -p synapse/EAID_000084/embeddings")
system("mkdir -p synapse/EAID_000085/embeddings")
system("mkdir -p synapse/EAID_000086/embeddings")
system("mkdir -p synapse/EAID_000087/embeddings")
system("mkdir -p synapse/EAID_000088/embeddings")
system("mkdir synapse/EAID_000084/markers")
system("mkdir synapse/EAID_000085/markers")
system("mkdir synapse/EAID_000086/markers")
system("mkdir synapse/EAID_000087/markers")
system("mkdir synapse/EAID_000088/markers")
system("mkdir synapse/EAID_000084/labels")
system("mkdir synapse/EAID_000085/labels")
system("mkdir synapse/EAID_000086/labels")
system("mkdir synapse/EAID_000087/labels")
system("mkdir synapse/EAID_000088/labels")
system("mkdir -p synapse/EAID_000084/figures/snrna")
system("mkdir -p synapse/EAID_000085/figures/snrna")
system("mkdir -p synapse/EAID_000086/figures/snrna")
system("mkdir -p synapse/EAID_000087/figures/snrna")
system("mkdir -p synapse/EAID_000088/figures/snrna")
system("mkdir -p synapse/EAID_000084/figures/snatac")
system("mkdir -p synapse/EAID_000085/figures/snatac")
system("mkdir -p synapse/EAID_000086/figures/snatac")
system("mkdir -p synapse/EAID_000087/figures/snatac")
system("mkdir -p synapse/EAID_000088/figures/snatac")
system("mkdir -p synapse/EAID_000084/auxiliary_data/snrna")
system("mkdir -p synapse/EAID_000085/auxiliary_data/snrna")
system("mkdir -p synapse/EAID_000086/auxiliary_data/snrna")
system("mkdir -p synapse/EAID_000087/auxiliary_data/snrna")
system("mkdir -p synapse/EAID_000088/auxiliary_data/snrna")
system("mkdir -p synapse/EAID_000084/auxiliary_data/snatac")
system("mkdir -p synapse/EAID_000085/auxiliary_data/snatac")
system("mkdir -p synapse/EAID_000086/auxiliary_data/snatac")
system("mkdir -p synapse/EAID_000087/auxiliary_data/snatac")
system("mkdir -p synapse/EAID_000088/auxiliary_data/snatac")

# Functions

## Reverse complement function

rev.comp<-function(x,rev=TRUE)
{
x<-toupper(x)
y<-rep("N",nchar(x))
xx<-unlist(strsplit(x,NULL))
for (bbb in 1:nchar(x))
    {
        if(xx[bbb]=="A") y[bbb]<-"T"        
        if(xx[bbb]=="C") y[bbb]<-"G"        
        if(xx[bbb]=="G") y[bbb]<-"C"        
        if(xx[bbb]=="T") y[bbb]<-"A"
    }
if(rev==FALSE) 
    {
    for(ccc in (1:nchar(x)))
        {
        if(ccc==1) yy<-y[ccc] else yy<-paste(yy,y[ccc],sep="")
        }
    }
if(rev==T)
    {
    zz<-rep(NA,nchar(x))
    for(ccc in (1:nchar(x)))
        {
        zz[ccc]<-y[nchar(x)+1-ccc]
        if(ccc==1) yy<-zz[ccc] else yy<-paste(yy,zz[ccc],sep="")
        }
    }
    return(yy)    
}

## Get tissue-level celltype metadata  

get_metadata = function(path){
    df = read.csv(path,row.names ="X")
    df = df[,c("cellID","timepoint","sex","tissue","sample","seurat_clusters",
               "gen_celltype","celltypes","subtypes")]
}

## Get raw ENCODE data matrices, make Seurat object, and pull metadata from object

get_orig_counts = function(file, expt_metadata){
    expt_metadata = expt_metadata[expt_metadata$file_accession == file,]
    
    if (unique(expt_metadata$technology)=="10x"){
        counts = readMM(paste0("snrna/counts_10x/",file,"/GeneFull_Ex50pAS/raw/UniqueAndMult-EM.mtx"))
        barcodes = read.delim(paste0("snrna/counts_10x/",file,"/GeneFull_Ex50pAS/raw/barcodes.tsv"),header = F, col.names="barcode")
        features = read.delim(paste0("snrna/counts_10x/",file,"/GeneFull_Ex50pAS/raw/features.tsv"),header = F) 
        colnames(counts) = colnames(counts) = paste0(barcodes$barcode,".",expt_metadata$library_accession)
        rownames(counts) = features$V2
        out = counts
        
        } else {
        counts = readMM(paste0("snrna/counts_parse/",file,"/GeneFull_Ex50pAS/raw/UniqueAndMult-EM.mtx"))
        barcodes = read.delim(paste0("snrna/counts_parse/",file,"/GeneFull_Ex50pAS/raw/barcodes.tsv"),header = F, col.names="barcode")
        features = read.delim(paste0("snrna/counts_parse/",file,"/GeneFull_Ex50pAS/raw/features.tsv"),header = F) 
        colnames(counts) = colnames(counts) = paste0(barcodes$barcode,".",expt_metadata$library_accession)
        rownames(counts) = features$V2
        out = counts
        
        }

}

make_seurat_obj = function(tissue, expt_metadata){
    expt_metadata = expt_metadata[expt_metadata$tissue == tissue,]
    files = expt_metadata$file_accession
    
    counts = list()
    for (j in 1:length(files)){
    counts[[j]] = get_orig_counts(files[j], expt_metadata)
    }
    counts_mat = do.call(cbind, counts)
    obj = CreateSeuratObject(counts = counts_mat, min.cells = 0, min.features = 0)
    obj[["percent.mt"]] = PercentageFeatureSet(obj, pattern = "^mt-")
    obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa")
    out = obj

}

format_rna_metadata = function(obj,expt_metadata,celltype_metadata){
    cells = obj@meta.data
    cells$cellID = rownames(cells)
    cells$rna_barcode= do.call("rbind", strsplit(as.character(cells$cellID), "[.]"))[,1]
    cells$library_accession= do.call("rbind", strsplit(as.character(cells$cellID), "[.]"))[,2]
    cells = left_join(cells, expt_metadata)
    cells = left_join(cells, celltype_metadata, by = "cellID")
    cells$passed_filtering[is.na(cells$passed_filtering)] = 0
    cells$nUMI = round(colSums(obj@assays$RNA@counts), digits = 0)
    cells = cells[,c("cellID","experiment_accession",
                                 "rna_barcode","atac_experiment_accession",
                                 "nUMI","percent.mt","percent.ribo","passed_filtering")]
    out = cells
}

## Get ATAC metadata from Pre-Filter-Metadata.rds files

get_archr_metadata = function(atac_experiment_metadata, tissue){
    atac_experiment_metadata = atac_experiment_metadata[atac_experiment_metadata$tissue == tissue,]
    libs = atac_experiment_metadata$rna_library_accession
    
    atac = list()
    for (i in 1:length(libs)){
        libID = libs[i]
        path = paste0("snatac/archr/QualityControl/",libID,"/",libID,"-Pre-Filter-Metadata.rds")
        atac[[i]] = readRDS(path)
        atac[[i]]$rna_library_accession = libID
    }
    cells = as.data.frame(do.call(rbind, atac))
    # read in 10x barcode sequences
    atac_bcs = read.delim("snatac/ref/atac_737K-arc-v1.txt",header = F)
    rna_bcs = read.delim("snatac/ref/gene_exp_737K-arc-v1.txt",header = F)
    bcs = cbind(atac_bcs,rna_bcs)
    colnames(bcs) = c("atac_bc","rna_bc")
    
    # match RNA and ATAC cell barcodes
    cells$original_atac_bc = do.call("rbind", strsplit(as.character(cells$cellNames), "#"))[,2]
    cells$atac_bc = sapply(cells$original_atac_bc,rev.comp)
    cells = dplyr::left_join(cells, bcs) # merge by atac_bc
    cells$cellID = paste0(cells$rna_bc,".",cells$rna_library_accession) # re-create RNA cell IDs
    out = cells
}


## Join ATAC and RNA metadata, save file

save_cell_metadata = function(rna_metadata, atac_metadata, eaid){
    full_metadata = full_join(rna_metadata, atac_metadata, by = "cellID")
    full_metadata = full_metadata[,c("cellID",
                                     "experiment_accession","rna_barcode",
                                     "atac_experiment_accession","original_atac_bc",
                                     "nUMI","nFrags",
                                     "percent.mt","percent.ribo","TSSEnrichment",
                                     "passed_filtering")]
    colnames(full_metadata) = c("cell_id","rna_dataset","rna_barcode",
                                "atac_dataset","atac_barcode",
                                "rna_umi_count","atac_fragment_count",
                                "rna_frac_mito","rna_frac_ribo","atac_tss_enrichment",
                                "passed_filtering")
    fn = paste0("synapse/",eaid,"/metadata_noheader.tsv")
    write.table(full_metadata,file=fn,
            sep="\t",
            quote=F,
            row.names = F)
    out = paste0("synapse/",eaid,"/metadata.tsv")
    system(paste0("cat synapse/header/cell_metadata_header ", fn, "> ", out))
    system(paste0("rm ", fn))
    system(paste0("gzip ", out))
}

## Get RNA embeddings

get_rna_embedding = function(rna_obj, eaid){
    rna_umap = data.frame(cell_id = rna_obj$cellID, 
                      UMAP_1 = rna_obj@reductions$umap@cell.embeddings[,1],
                      UMAP_2 = rna_obj@reductions$umap@cell.embeddings[,2])
    fn = paste0("synapse/",eaid,"/embeddings/rna_umap_noheader.tsv")
    write.table(rna_umap,file=fn,
            sep="\t",
            quote=F,
            row.names = F)
    out = paste0("synapse/",eaid,"/embeddings/rna_umap.tsv")
    system(paste0("cat synapse/header/rna_embedding_coordinates_header ", fn, "> ", out))
    system(paste0("rm ", fn))
    system(paste0("gzip ", out))
}

## Get ATAC embeddings

get_atac_embedding = function(archr_obj, eaid){
    coords = getEmbedding(archr_obj)
    atac_umap = data.frame(cell_id = archr_obj$cellID, 
                      UMAP_1 = coords[,1],
                      UMAP_2 = coords[,2])
    
    fn = paste0("synapse/",eaid,"/embeddings/atac_umap_noheader.tsv")
    write.table(atac_umap,file=fn,
            sep="\t",
            quote=F,
            row.names = F)
    out = paste0("synapse/",eaid,"/embeddings/atac_umap.tsv")
    system(paste0("cat synapse/header/atac_embedding_coordinates_header ", fn, "> ", out))
    system(paste0("rm ", fn))
    system(paste0("gzip ", out))
}


## Get marker genes

## Call cell type marker genes 
get_markers = function(combined.sct,eaid){
    geneinfo = read.csv("snrna/ref/gene_id_to_name.csv")
    geneinfo$gene = make.unique(geneinfo$gene_name)
    
    DefaultAssay(combined.sct)= "RNA"
    Idents(combined.sct) = "celltypes"
    markers <- FindAllMarkers(combined.sct, 
                             only.pos = TRUE, 
                             min.pct = 0.25, 
                             logfc.threshold = 0.25, 
                             verbose = T)
    markers=left_join(markers,geneinfo)
    markers$is_enriched = 1
    markers = markers[,c("gene_id","gene_name","is_enriched","avg_log2FC","p_val_adj","pct.1","pct.2","cluster")]
    markers_list <- split(markers, f = markers$cluster) 
    
    for (i in 1:length(markers_list)){
        celltype = unique(markers_list[[i]]$cluster)
        markers_list[[i]]$cluster = NULL
        
        celltype = gsub("[/]",".",celltype)
        fn = paste0("synapse/",eaid,"/markers/",celltype,"_markersnoheader.tsv")
        write.table(markers_list[[i]],file=fn,
            sep="\t",
            quote=F,
            row.names = F)
        out = paste0("synapse/",eaid,"/markers/",celltype,"_markers.tsv")
        system(paste0("cat synapse/header/markers_header ", fn, "> ", out))
        system(paste0("rm ", fn))
        system(paste0("gzip ", out))
        }

}

## Make cell type labels file

cell_type_labels = function(rna_obj, eaid){
    celltype_meta = rna_obj@meta.data
    celltype_meta = celltype_meta[,c("cellID","celltypes")]
    celltype_meta$cell_type_id = celltype_meta$celltypes
    
    for (i in 1:length(unique(celltype_meta$celltypes))){
        celltype = unique(celltype_meta$celltypes)[i]
        celltype = gsub("[/]",".",celltype)
        fn = paste0("synapse/",eaid,"/markers/",celltype,"_markers.tsv.gz")
        markers = fread(fn)
        name = paste0(celltype, "_membership_score")
        rna_obj[[name]] = PercentageFeatureSet(rna_obj, features = make.unique(markers$gene_name),  assay ="RNA")
        
        celltype = unique(celltype_meta$celltypes)[i]
        rna_celltype_obj = subset(rna_obj, subset = celltypes == celltype)
        celltype_meta$membership_score[celltype_meta$celltypes == celltype] = rna_celltype_obj[[name]][,1]
   }
    celltype_meta$cell_type_id = celltype_meta$celltypes
    colnames(celltype_meta) = c("cell_id","cell_type_id","cell_type_name","membership_score")
    rownames(celltype_meta) = NULL
    
    fn = paste0("synapse/",eaid,"/labels/cell_type_labels_noheader.tsv")
    write.table(celltype_meta,file=fn,
                sep="\t",
                quote=F,
                row.names = F)
    out = paste0("synapse/",eaid,"/labels/cell_type_labels.tsv")
    system(paste0("cat synapse/header/labels_header ", fn, "> ", out))
    system(paste0("rm ", fn))
    system(paste0("gzip ", out))
}

## Collect figures

collect_figures = function(tissue, eaid, archr_proj_name){
    setwd("/share/crsp/lab/seyedam/share/enc4_mouse/")
    system(paste0("cp -r snrna/plots/",tissue,"/* synapse/", eaid,"/figures/snrna/"))
    system(paste0("cp snatac/archr/",archr_proj_name,"/Plots/* synapse/", eaid,"/figures/snatac"))
    setwd(paste0("synapse/", eaid))
    system(paste0("tar -czf figures.tar.gz figures/ --remove-files"))
    setwd("/share/crsp/lab/seyedam/share/enc4_mouse/")
}

## Collect data

collect_data = function(tissue, eaid, archr_proj_name){
    setwd("/share/crsp/lab/seyedam/share/enc4_mouse/")
    setwd(paste0("synapse/", eaid))
    system(paste0("cp ../../snrna/seurat/",tissue,"* auxiliary_data/snrna/"))
    system(paste0("cp ../../snrna/seurat/markers/",tissue,"_cluster* auxiliary_data/snrna/"))
    system(paste0("cp -r ../../snatac/archr/",archr_proj_name," auxiliary_data/snatac"))
    system(paste0("tar -czf auxiliary_data.tar.gz auxiliary_data/ --remove-files"))
    setwd("/share/crsp/lab/seyedam/share/enc4_mouse/")
}

# Make datasets.txt files

expt_metadata = read.delim("snrna/ref/enc4_mouse_snrna_metadata.tsv")
atac_experiment_metadata = read.delim("snatac/ref/enc4_mouse_snatac_metadata.tsv")

datasets = as.data.frame(c(expt_metadata[expt_metadata$tissue == "Adrenal","experiment_accession"],
             atac_experiment_metadata[atac_experiment_metadata$tissue == "Adrenal","experiment_accession"]))
colnames(datasets) = "datasets"
write.table(datasets, file="synapse/EAID_000084/datasets.txt",quote=F,row.names = F, col.names = F)                

datasets = as.data.frame(c(expt_metadata[expt_metadata$tissue == "Cortex","experiment_accession"],
             atac_experiment_metadata[atac_experiment_metadata$tissue == "Cortex","experiment_accession"]))
colnames(datasets) = "datasets"
write.table(datasets, file="synapse/EAID_000085/datasets.txt",quote=F,row.names = F, col.names = F)                

datasets = as.data.frame(c(expt_metadata[expt_metadata$tissue == "Hippocampus","experiment_accession"],
             atac_experiment_metadata[atac_experiment_metadata$tissue == "Hippocampus","experiment_accession"]))
colnames(datasets) = "datasets"
write.table(datasets, file="synapse/EAID_000086/datasets.txt",quote=F,row.names = F, col.names = F)                

datasets = as.data.frame(c(expt_metadata[expt_metadata$tissue == "Heart","experiment_accession"],
             atac_experiment_metadata[atac_experiment_metadata$tissue == "Heart","experiment_accession"]))
colnames(datasets) = "datasets"
write.table(datasets, file="synapse/EAID_000087/datasets.txt",quote=F,row.names = F, col.names = F)                

datasets = as.data.frame(c(expt_metadata[expt_metadata$tissue == "Gastrocnemius","experiment_accession"],
             atac_experiment_metadata[atac_experiment_metadata$tissue == "Gastrocnemius","experiment_accession"]))
colnames(datasets) = "datasets"
write.table(datasets, file="synapse/EAID_000088/datasets.txt",quote=F,row.names = F, col.names = F)                

# Cell Metadata File

expt_metadata = read.delim("snrna/ref/enc4_mouse_snrna_metadata.tsv")
atac_experiment_metadata = read.delim("snatac/ref/enc4_mouse_snatac_metadata.tsv")

## Adrenal gland
print("Saving adrenal metadata...")
obj = make_seurat_obj("Adrenal",expt_metadata)
rna_metadata = format_rna_metadata(obj,expt_metadata, celltype_metadata)

print(table(rna_metadata$passed_filtering))

atac_metadata = get_archr_metadata(atac_experiment_metadata, "Adrenal")

save_cell_metadata(rna_metadata, atac_metadata, "EAID_000084")



## Left cortex
print("Saving cortex metadata...")
obj = make_seurat_obj("Cortex",expt_metadata)
rna_metadata = format_rna_metadata(obj,expt_metadata, celltype_metadata)

print(table(rna_metadata$passed_filtering))

atac_metadata = get_archr_metadata(atac_experiment_metadata, "Cortex")

save_cell_metadata(rna_metadata, atac_metadata, "EAID_000085")


## Hippocampus
print("Saving hippocampus metadata...")
obj = make_seurat_obj("Hippocampus",expt_metadata)
rna_metadata = format_rna_metadata(obj,expt_metadata, celltype_metadata)

print(table(rna_metadata$passed_filtering))

atac_metadata = get_archr_metadata(atac_experiment_metadata, "Hippocampus")

save_cell_metadata(rna_metadata, atac_metadata, "EAID_000086")

## Heart
print("Saving heart metadata...")
obj = make_seurat_obj("Heart",expt_metadata)
rna_metadata = format_rna_metadata(obj,expt_metadata, celltype_metadata)

print(table(rna_metadata$passed_filtering))

atac_metadata = get_archr_metadata(atac_experiment_metadata, "Heart")

save_cell_metadata(rna_metadata, atac_metadata, "EAID_000087")

## Gastrocnemius
print("Saving gastroc metadata...")
obj = make_seurat_obj("Gastrocnemius",expt_metadata)
rna_metadata = format_rna_metadata(obj,expt_metadata, celltype_metadata)

print(table(rna_metadata$passed_filtering))

atac_metadata = get_archr_metadata(atac_experiment_metadata, "Gastrocnemius")

save_cell_metadata(rna_metadata, atac_metadata, "EAID_000088")

# Embeddings

## Adrenal gland

adr_rna = readRDS("snrna/seurat/adrenal_Parse_10x_integrated.rds")
get_rna_embedding(adr_rna, "EAID_000084")
adr_atac = loadArchRProject("snatac/archr/ENC4_Mouse_Adrenal/", showLogo = F)
get_atac_embedding(adr_atac, "EAID_000084")

## Left cortex

ctx_rna = readRDS("snrna/seurat/cortex_Parse_10x_integrated.rds")
get_rna_embedding(ctx_rna, "EAID_000085")
ctx_atac = loadArchRProject("snatac/archr/ENC4_Mouse_Cortex/", showLogo = F)
get_atac_embedding(ctx_atac, "EAID_000085")

## Hippocampus

hc_rna = readRDS("snrna/seurat/hippocampus_Parse_10x_integrated.rds")
get_rna_embedding(hc_rna, "EAID_000086")
hc_atac = loadArchRProject("snatac/archr/ENC4_Mouse_Hippocampus/", showLogo = F)
get_atac_embedding(hc_atac, "EAID_000086")

## Heart

hrt_rna = readRDS("snrna/seurat/heart_Parse_10x_integrated.rds")
get_rna_embedding(hrt_rna, "EAID_000087")
hrt_atac = loadArchRProject("snatac/archr/ENC4_Mouse_Heart/", showLogo = F)
get_atac_embedding(hrt_atac, "EAID_000087")

## Gastrocnemius

gas_rna = readRDS("snrna/seurat/gastrocnemius_Parse_10x_integrated.rds")
get_rna_embedding(gas_rna, "EAID_000088")
gas_atac = loadArchRProject("snatac/archr/ENC4_Mouse_Gastrocnemius/", showLogo = F)
get_atac_embedding(gas_atac, "EAID_000088")

# Marker genes
get_markers(adr_rna,"EAID_000084")
get_markers(ctx_rna,"EAID_000085")
get_markers(hc_rna,"EAID_000086")
get_markers(hrt_rna,"EAID_000087")
get_markers(gas_rna,"EAID_000088")

# Cell Type Labels
cell_type_labels(adr_rna, "EAID_000084")
cell_type_labels(ctx_rna, "EAID_000085")
cell_type_labels(hc_rna, "EAID_000086")
cell_type_labels(hrt_rna, "EAID_000087")
cell_type_labels(gas_rna, "EAID_000088")

# Collect figures into tar zipped folder
collect_figures("adrenal","EAID_000084","ENC4_Mouse_Adrenal")
collect_figures("cortex","EAID_000085","ENC4_Mouse_Cortex")
collect_figures("hippocampus","EAID_000086","ENC4_Mouse_Hippocampus")
collect_figures("heart","EAID_000087","ENC4_Mouse_Heart")
collect_figures("gastrocnemius","EAID_000088","ENC4_Mouse_Gastrocnemius")

# Collect auxiliary data into tar zipped folder
# Including Seurat object, Seurat object metadata as a csv file, cluster-level markers genes tsv file, and ArchR project folder.
collect_data("adrenal","EAID_000084","ENC4_Mouse_Adrenal")
collect_data("cortex","EAID_000085","ENC4_Mouse_Cortex")
collect_data("hippocampus","EAID_000086","ENC4_Mouse_Hippocampus")
collect_data("heart","EAID_000087","ENC4_Mouse_Heart")
collect_data("gastrocnemius","EAID_000088","ENC4_Mouse_Gastrocnemius")
