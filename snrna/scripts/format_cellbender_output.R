library(Matrix)
library(Seurat)

setwd("/share/crsp/lab/seyedam/share/enc4_mouse/snrna/")
metadata = read.delim("ref/enc4_mouse_snrna_metadata.tsv")
metadata = metadata[metadata$technology == "10x",]

files = metadata$file_accession
batches = unique(metadata$experiment_batch)

for (i in 1:length(batches)){
    metadata_expt = metadata[metadata$experiment_batch == batches[i],]
    files = metadata_expt$file_accession

    counts = Read10X_h5(paste0("counts_10x/",files,"/cellbender.h5"),unique.features = FALSE)
    bcs = read.delim(paste0("counts_10x/",files,"/barcodes.tsv"),header=F)
    colnames(counts) = bcs$V1
    
    counts = counts[,colSums(counts) > 500] # > 500 UMI

    writeMM(counts,file=paste0("scrublet/",batches[i],"_matrix.mtx"))
    write.table(colnames(counts),file=paste0("scrublet/",batches[i],"_barcodes.tsv"),quote=F,row.names=F,col.names=F,sep="\t")
    write.table(rownames(counts),file=paste0("scrublet/",batches[i],"_genes.tsv"),quote=F,row.names=F,col.names=F,sep="\t")
}

