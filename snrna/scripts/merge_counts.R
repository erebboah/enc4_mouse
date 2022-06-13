# Load and format data
library(Matrix)
library(stringr)
library(Seurat)

#Read in gene reference and filter by biotypes of interst
setwd("/share/crsp/lab/seyedam/share/enc4_mouse/snrna/")
metadata = read.delim("ref/enc4_mouse_snrna_metadata.tsv")

gene_id_to_name = read.csv("ref/gene_id_to_name.csv")

gene_id_to_name = gene_id_to_name[gene_id_to_name$gene_type %in% c("protein_coding",
                                                                      "miRNA",
                                                                      # lncRNA https://www.gencodegenes.org/pages/biotypes.html
                                                                      "3prime_overlapping_ncRNA",
                                                                      "antisense",
                                                                      "bidirectional_promoter_lncRNA",
                                                                      "lincRNA",
                                                                      "macro_lncRNA",
                                                                      "non_coding",
                                                                      "processed_transcript",
                                                                      "sense_intronic",
                                                                      "sense_overlapping",
                                                                      # pseudogenes
                                                                      "processed_pseudogene",
                                                                      "transcribed_unprocessed_pseudogene",
                                                                      "unprocessed_pseudogene",
                                                                      "transcribed_processed_pseudogene",
                                                                      "unitary_pseudogene",
                                                                      "pseudogene",
                                                                      "polymorphic_pseudogene",
                                                                      "transcribed_unitary_pseudogene",
                                                                      "translated_unprocessed_pseudogene"),]


# Read in data and filter
## 1. gene with biotype specified above 
## 2. > 500 UMIs/nucleus

get_counts = function(file){
    metadata = metadata[metadata$file_accession == file,]
    path = paste0("counts/",file,"/GeneFull_Ex50pAS/raw/")
    barcodes = read.delim(paste0(path,'barcodes.tsv'),header = F, col.names="barcode")
    counts = readMM(paste0(path,'UniqueAndMult-EM.mtx'))
    features = read.delim(paste0(path,'features.tsv'),header = F) 
    rownames(counts) = features$V2 
    colnames(counts) = paste0(barcodes$barcode,".",metadata$library_accession) # append library accession to cell barcode
    counts = counts[rownames(counts) %in% gene_id_to_name$gene_name,]
    counts = counts[,colSums(counts) > 500] # > 500 UMI
    out = counts
}


#make seurat objects
seurat_obj = function(counts,metadata){
    obj = CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0)
    obj@meta.data = cbind(obj@meta.data,metadata)
    obj[["percent.mt"]] = PercentageFeatureSet(obj, pattern = "^mt-")
    obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa")
    out = obj
}

tissues = c("Adrenal","Gastrocnemius","Heart","Cortex","Hippocampus")

for (k in 1:length(tissues)){
    tissue = tissues[k]
    files = metadata$file_accession[metadata$tissue == tissue]

    matrix_list = list()
    for (i in 1:length(files)){
        matrix_list[[i]] = get_counts(files[i])
    }

    matrix = matrix_list[[1]]
    for (j in 2:length(files)){
        matrix = RowMergeSparseMatrices(matrix,matrix_list[[j]])
    }

    tissue_metadata = data.frame(cellID = colnames(matrix))
    tissue_metadata$library_accession = do.call("rbind", strsplit(tissue_metadata$cellID, "[.]"))[,2]
    tissue_metadata=dplyr::left_join(tissue_metadata, metadata, by = "library_accession")

    writeMM(matrix,file=paste0("counts/",tissue,"_counts_500UMI.mtx"))
    write.csv(tissue_metadata,file=paste0("counts/",tissue,"_counts_500UMI_metadata.csv"))
    write.csv(rownames(matrix),file=paste0("counts/",tissue,"_counts_500UMI_genes.csv"))
}

filters = list()
filters[[1]] = c()

tissue_obj_list = list()
for (i in 1:length(tissues)){
    tissue = tissues[i]
    counts = readMM(paste0("counts/",tissue,"_counts_500UMI.mtx"))
    genes = read.csv(paste0("counts/",tissue,"_counts_500UMI_genes.csv"))
    metadata = read.csv(paste0("counts/",tissue,"_counts_500UMI_metadata.csv"))
    genes = read.csv(paste0("counts/",tissue,"_counts_500UMI_genes.csv"),row.names=1)
    colnames(counts) = metadata$cellID
    rownames(counts) = genes[,1]
    tissue_obj_list[[i]] = seurat_obj(counts,metadata)
    saveRDS(tissue_obj_list[[i]], file = paste0("seurat/",tissue,"_seurat.rds"))
    
    # split for scrublet processing

    for (k in 1:length(tissue_obj_list)){
    obj_list = list()
    obj_list = SplitObject(tissue_obj_list[[k]], split.by = "experiment_batch")

        for (j in 1:length(obj_list)){
        counts = obj_list[[j]]@assays$RNA@counts
        meta  = obj_list[[j]]@meta.data
        writeMM(counts,file=paste0("scrublet/",tissue,"_",unique(meta$experiment_batch),"_counts_500UMI.mtx"))
        write.csv(meta,file=paste0("scrublet/",tissue,"_",unique(meta$experiment_batch),"_counts_500UMI_metadata.csv"))
        }
    }

}

