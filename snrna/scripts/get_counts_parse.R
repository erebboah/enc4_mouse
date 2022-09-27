# Load and format data
library(Matrix)
library(stringr)
library(Seurat)

#Read in gene reference and filter by biotypes of interst
setwd("../../snrna/")
metadata = read.delim("ref/enc4_mouse_snrna_metadata.tsv")
metadata = metadata[metadata$technology == "Parse",]

gene_id_to_name = read.delim("ref/geneInfo.tab",header=F,col.names=c("gene_id","gene_name","gene_type"))
gene_id_to_name = gene_id_to_name[-1,]

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
    
    if (file %in% c("ENCFF907ZKD","ENCFF006DTQ")) {
    path = paste0("counts_parse/",file,"/Gene/raw/")
    barcodes = read.delim(paste0(path,'barcodes.tsv'),header = F, col.names="barcode")
    counts = readMM(paste0(path,'UniqueAndMult-EM.mtx'))
    features = read.delim(paste0(path,'features.tsv'),header = F) 
    rownames(counts) = features$V2 
    colnames(counts) = paste0(barcodes$barcode,".",metadata$library_accession) # append library accession to cell barcode
    counts = counts[,colSums(counts) > 500] # > 500 UMI
    counts = counts[rownames(counts) %in% gene_id_to_name$gene_name,]
    out = counts    
    } else {
    
    path = paste0("counts_parse/",file,"/GeneFull_Ex50pAS/raw/")
    barcodes = read.delim(paste0(path,'barcodes.tsv'),header = F, col.names="barcode")
    counts = readMM(paste0(path,'UniqueAndMult-EM.mtx'))
    features = read.delim(paste0(path,'features.tsv'),header = F) 
    rownames(counts) = features$V2 
    colnames(counts) = paste0(barcodes$barcode,".",metadata$library_accession) # append library accession to cell barcode
    counts = counts[,colSums(counts) > 500] # > 500 UMI
    counts = counts[rownames(counts) %in% gene_id_to_name$gene_name,]
    out = counts
    }
}

files = metadata$file_accession
batches = unique(metadata$experiment_batch)


for (i in 1:length(batches)){
	metadata_expt = metadata[metadata$experiment_batch == batches[i],]
	files = metadata_expt$file_accession

	counts = get_counts(files[1])

	for (j in 2:length(files)){
	counts_adding = get_counts(files[j])
	counts = cbind(counts, counts_adding)
}
        writeMM(counts,file=paste0("scrublet/",batches[i],"_matrix.mtx"))
        write.table(colnames(counts),file=paste0("scrublet/",batches[i],"_barcodes.tsv"),quote=F,row.names=F,col.names=F,sep="\t")
        write.table(rownames(counts),file=paste0("scrublet/",batches[i],"_genes.tsv"),quote=F,row.names=F,col.names=F,sep="\t")
}

gene_id_to_name_save = gene_id_to_name[gene_id_to_name$gene_name %in% rownames(counts),]
gene_id_to_name_save = gene_id_to_name_save[match(rownames(counts), gene_id_to_name_save$gene_name),]
write.csv(gene_id_to_name_save,file="ref/gene_id_to_name.csv")



