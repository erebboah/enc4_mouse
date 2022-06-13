library(Matrix)
library(stringr)
library(Seurat)

setwd("/share/crsp/lab/seyedam/share/enc4_mouse/snrna/")
metadata = read.delim("ref/enc4_mouse_snrna_metadata.tsv")

scrub_results = dir("scrublet")
scrub_results = scrub_results[grep(".csv",scrub_results)]
new_metadata = read.csv(paste0("scrublet/",scrub_results[1]))
for (i in 2:length(scrub_results)){
    new_metadata_adding = read.csv(paste0("scrublet/",scrub_results[i]))
    new_metadata = rbind(new_metadata, new_metadata_adding)
}

all_filters = list()
# lower UMI, upper UMI, lower gene, lower doublet, lower mito
# adrenal
all_filters[[1]] = list()
all_filters[[1]][[1]]= c(1000,30000,500,0.25,1) 
all_filters[[1]][[2]] = c(500,30000,500,0.25,1) 
all_filters[[1]][[3]] = c(1000,50000,500,0.25,1) 

# gastroc
all_filters[[2]] = list()
all_filters[[2]][[1]]= c(500,30000,300,0.25,5) 
all_filters[[2]][[2]] = c(500,30000,300,0.25,5) 
all_filters[[2]][[3]] = c(1000,50000,300,0.25,5) 

# heart
all_filters[[3]] = list()
all_filters[[3]][[1]]= c(1000,30000,500,0.25,1) 
all_filters[[3]][[2]] = c(500,30000,500,0.25,1) 
all_filters[[3]][[3]] = c(1000,50000,500,0.25,1) 

# cortex
all_filters[[4]] = list()
all_filters[[4]][[1]]= c(1000,50000,500,0.2,0.5) 
all_filters[[4]][[2]] = c(500,30000,500,0.2,0.5) 
all_filters[[4]][[3]] = c(1000,50000,500,0.2,0.5) 

# HC
all_filters[[5]] = list()
all_filters[[5]][[1]]= c(1000,50000,500,0.2,0.5) 
all_filters[[5]][[2]] = c(500,30000,500,0.2,0.5) 
all_filters[[5]][[3]] = c(1000,50000,500,0.2,0.5) 


tissues = c("Adrenal","Gastrocnemius","Heart","Cortex","Hippocampus")

for (j in 1:length(tissues)){
    
tissue = tissues[j]

filters = all_filters[[j]]
    
obj_list = list()
obj = readRDS(paste0("seurat/",tissue,"_seurat.rds"))
new_metadata_obj = new_metadata[new_metadata$cellID %in% colnames(obj),]
new_metadata_obj = new_metadata_obj[match(colnames(obj), new_metadata_obj$cellID),]
obj$doublet_scores = new_metadata_obj$doublet_scores
obj$doublets = new_metadata_obj$doublets
# split by technology_depth...3 objects: 10x, parse shallow, and parse deep.
obj$technology_depth = paste0(obj$technology,"_",obj$depth1)
obj_list = SplitObject(obj, split.by = "technology_depth")

    for (i in 1:length(filters)){
    obj_list[[i]] <- subset(obj_list[[i]], 
                        subset = nCount_RNA > filters[[i]][1] & 
                        nCount_RNA < filters[[i]][2] & 
                        nFeature_RNA > filters[[i]][3] & 
                        doublet_scores < filters[[i]][4] & 
                        percent.mt < filters[[i]][5])
    }
   
obj_list <- lapply(X = obj_list, FUN = SCTransform, vars.to.regress = c("percent.mt","nFeature_RNA"))
features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 3000)
obj_list <- PrepSCTIntegration(object.list = obj_list, anchor.features = features)  
names(obj_list) = c("10x","shallow","deep")
# use Parse shallow as reference dataset
reference_dataset <- which(names(obj_list) == "shallow")
anchors <- FindIntegrationAnchors(object.list = obj_list, normalization.method = "SCT", anchor.features = features, reference = reference_dataset)
combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
combined.sct <- RunPCA(combined.sct, verbose = T, npcs = 50)
combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:30,spread=2)
combined.sct <- FindNeighbors(combined.sct, reduction = "pca", dims = 1:30)
combined.sct <- FindClusters(combined.sct)            
saveRDS(combined.sct,file=paste0(paste0("seurat/",tissue,"_Parse_10x_integrated.rds")))    
} 

