# Load and format data
library(Matrix)
library(stringr)
library(Seurat)

#Read in gene reference and filter by biotypes of interst
setwd("../../snrna/")
metadata = read.delim("ref/enc4_mouse_snrna_metadata.tsv")
metadata = metadata[metadata$technology == "10x",]

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
## 2. > 0 UMIs/nucleus (remove nuclei with 0 counts)

get_counts = function(file){
    metadata = metadata[metadata$file_accession == file,]
    
    path = paste0("counts_10x/",file,"/GeneFull_Ex50pAS/raw/")
    barcodes = read.delim(paste0(path,'barcodes.tsv'),header = F, col.names="barcode")
    counts = readMM(paste0(path,'UniqueAndMult-EM.mtx'))
    features = read.delim(paste0(path,'features.tsv'),header = F) 
    rownames(counts) = features$V2 
    colnames(counts) = paste0(barcodes$barcode,".",metadata$library_accession) # append library accession to cell barcode
    counts = counts[rownames(counts) %in% gene_id_to_name$gene_name,]
    counts = counts[,colSums(counts) > 0] # > 0 UMI
    out = counts
}

files = metadata$file_accession

for (i in 1:length(files)){
        counts = get_counts(files[i])

	writeMM(counts,file=paste0("counts_10x/",files[i],"/matrix.mtx"))
        write.table(colnames(counts),file=paste0("counts_10x/",files[i],"/barcodes.tsv"),quote=F,row.names=F,col.names=F,sep="\t")
        write.table(rownames(counts),file=paste0("counts_10x/",files[i],"/genes.tsv"),quote=F,row.names=F,col.names=F,sep="\t")
}
