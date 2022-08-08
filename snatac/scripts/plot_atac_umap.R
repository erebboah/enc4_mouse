suppressPackageStartupMessages(library(ArchR))
library(parallel)
library(viridis)
suppressPackageStartupMessages(library(Seurat))
set.seed(1234)
library(reshape2)
addArchRThreads(threads = 8) 
addArchRGenome("mm10")

plot_umap = function(atac_obj,
                   slot){
   
    
    ref = read.csv("/share/crsp/lab/seyedam/share/Heidi_Liz/ref/enc4_mouse_snrna_celltypes.csv",stringsAsFactors = F)
    ref$subtypes_color = ref$subtype_color
    ref$celltypes_color = ref$celltype_color
    this_color = paste0(slot,"_color")
    
    st = sort(unique(meta[,slot]))
    ref = ref[ref[,slot] %in% st,]
    ref = ref[match(st,ref[,slot]),]
    ref_plot = unique(ref[,c(slot,this_color)])
    
    meta = getCellColData(proj) 
    #plot
    coords=as.data.frame(getEmbedding(atac_obj, embedding = "UMAP")) # grab UMAP coordinates of cells
    df <- data.frame(UMAP_1 = coords[,1],UMAP_2 = coords[,2])

    # add expression and enrichmen info
    df$viz = meta[,slot]
    
    p1 = ggplot(df, aes(UMAP_1, UMAP_2, colour = viz)) +
    geom_point(size=0.1) + 
    theme_classic() + scale_color_manual(values = ref_plot[,2]) +
    labs(title = slot) +
    theme(text = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)) + 
    theme(
    title = element_text(size = 10),
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank())+ 
    labs(colour = slot) + guides(colour = guide_legend(override.aes = list(size=10)))
    
    return(p1)
}


suppressPackageStartupMessages(library(ArchR))
library(parallel)
library(viridis)
suppressPackageStartupMessages(library(Seurat))
set.seed(1234)
library(reshape2)
addArchRThreads(threads = 8) 
addArchRGenome("mm10")

plot_umap = function(atac_obj,
                   slot){
   
    
    ref = read.csv("/share/crsp/lab/seyedam/share/Heidi_Liz/ref/enc4_mouse_snrna_celltypes.csv",stringsAsFactors = F)
    ref$subtypes_color = ref$subtype_color
    ref$celltypes_color = ref$celltype_color
    this_color = paste0(slot,"_color")
    
    st = sort(unique(meta[,slot]))
    ref = ref[ref[,slot] %in% st,]
    ref = ref[match(st,ref[,slot]),]
    ref_plot = unique(ref[,c(slot,this_color)])
    
    meta = getCellColData(proj) 
    #plot
    coords=as.data.frame(getEmbedding(atac_obj, embedding = "UMAP")) # grab UMAP coordinates of cells
    df <- data.frame(UMAP_1 = coords[,1],UMAP_2 = coords[,2])

    # add expression and enrichmen info
    df$viz = meta[,slot]
    
    p1 = ggplot(df, aes(UMAP_1, UMAP_2, colour = viz)) +
    geom_point(size=0.1) + 
    theme_classic() + scale_color_manual(values = ref_plot[,2]) +
    labs(title = slot) +
    theme(text = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)) + 
    theme(
    title = element_text(size = 10),
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank())+ 
    labs(colour = slot) + guides(colour = guide_legend(override.aes = list(size=10)))
    
    return(p1)
}


setwd("/share/crsp/lab/seyedam/share/enc4_mouse/snatac/archr/")

proj = loadArchRProject(path = "ENC4_Mouse/")

pdf(file="/share/crsp/lab/seyedam/share/enc4_mouse/snrna/archr/ENC4_Mouse/atac_umaps_correct_colors.pdf",width=12,height=12)
print(plot_umap(proj,"gen_celltype"))
print(plot_umap(proj,"celltypes"))
dev.off()

