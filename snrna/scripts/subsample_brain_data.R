library(Matrix)
library(Seurat)
library(dplyr)
set.seed("1234")

# subsample metadata
metadata = read.csv("metadata.csv")
celltypes = unique(metadata$cell_type_alias_label)

atlas_celltypes = as.data.frame(table(metadata$cell_type_alias_label))

atlas_celltypes_under1000 = atlas_celltypes$Var1[atlas_celltypes$Freq < 1000]
subsample_under1000 =  metadata[metadata$cell_type_alias_label %in% atlas_celltypes_under1000,]
atlas_celltypes_over1000 = atlas_celltypes$Var1[atlas_celltypes$Freq >= 1000]
subsample_over1000 =  metadata[metadata$cell_type_alias_label %in% atlas_celltypes_over1000,]
subsample_over1000 =  subsample_over1000 %>% group_by(cell_type_alias_label) %>% slice_sample(n=1000)
subsample = rbind(subsample_under1000,subsample_over1000)

write.csv(subsample,file="brain_atlas_metadata_1000_per_subtype_250734.csv",quote=F,row.names=F)

# load CSV files
for (i in 1:12){
  brain = read.csv(paste0("matrix_",i,".csv"),header=F)
  genes = as.character(as.vector(brain[1,2:ncol(brain)]))
  
  print(table(brain[,1] %in% subsample$sample_name))
  subsample_section = subsample[subsample$sample_name %in% brain[,1],]
  brain_subsample = brain[brain[,1] %in% subsample_section$sample_name,]
  subsample_section = subsample_section[match(brain_subsample[,1],subsample_section$sample_name),]
  
  # format as sparse matrix
  counts = brain_subsample 
  counts[,1] = NULL
  counts = t(counts)
  counts = Matrix(counts,sparse=T)
  
  colnames(counts) = subsample_section$sample_name
  rownames(counts) = genes
  
  writeMM(counts,file=paste0("matrix_",i,"_subsampled.mtx"))
  write.table(rownames(counts),file=paste0("matrix_",i,"_genes.tsv"),sep="\t",quote=F,row.names=F,col.names = F)
  write.table(colnames(counts),file=paste0("matrix_",i,"_barcodes.tsv"),sep="\t",quote=F,row.names=F,col.names = F)
}

# combine subsampled matrices for 1 sparse matrix
counts = list()
for (i in 1:12){
  counts[[i]] = readMM(paste0("matrix_",i,"_subsampled.mtx"))
  genes = read.delim(paste0("matrix_",i,"_genes.tsv"),header=F)
  barcodes = read.delim(paste0("matrix_",i,"_barcodes.tsv"),header=F)
  rownames(counts[[i]]) = genes$V1
  colnames(counts[[i]]) = barcodes$V1
}
brain_sub = RowMergeSparseMatrices(counts[[1]], c(counts[[2]], counts[[3]], counts[[4]], counts[[5]], 
                                                  counts[[6]], counts[[7]], counts[[8]], counts[[9]], counts[[10]],
                                                  counts[[11]], counts[[12]]))

writeMM(brain_sub,file="brain_atlas_metadata_1000_per_subtype_250734.mtx")
write.table(rownames(brain_sub),file=paste0("brain_atlas_metadata_1000_per_subtype_250734_genes.tsv"),
            sep="\t",quote=F,row.names=F,col.names = F)
write.table(colnames(brain_sub),file=paste0("brain_atlas_metadata_1000_per_subtype_250734_barcodes.tsv"),
            sep="\t",quote=F,row.names=F,col.names = F)

