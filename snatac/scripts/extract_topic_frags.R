suppressPackageStartupMessages(library(ArchR))
library(parallel)
library(viridis)
suppressPackageStartupMessages(library(Seurat))
set.seed(1234)
library(reshape2)
library(Repitools)
library(tidyverse)
addArchRThreads(threads = 32)
addArchRGenome("mm10")

setwd("/share/crsp/lab/seyedam/share/enc4_mouse/snatac/archr_filt/")

proj = loadArchRProject(path = "ENC4_Mouse_filt/")

topic_cells = read.csv("/share/crsp/lab/seyedam/share/enc4_mouse/snatac/ref/selectedcells_top5%_130cells.csv",row.names=1)
meta = getCellColData(proj)

atac_cells_filt = proj$cellNames[meta$tissue == "Heart"]

proj_heart = subsetArchRProject(
    ArchRProj = proj,
    outputDirectory = "ENC4_Mouse_heart_filt",
    cells = atac_cells_filt,
    dropCells = TRUE, #  drop cells that are not in ArchRProject from corresponding Arrow Files
    force = TRUE)

meta = getCellColData(proj_heart)

for (i in 1:ncol(topic_cells)) {
    this_topic_cells = unique(topic_cells[,i][topic_cells[,i] != ""])
    
    meta_subset = meta[meta$cellID %in% this_topic_cells,]
            
    topic_frags = getFragmentsFromProject(ArchRProj = proj_heart,
                            cellNames = rownames(meta_subset))
    df = list()
    for (j in 1:length(unique(proj_heart$Sample)){ # n arrow files
        df[[j]] <- annoGR2DF(topic_frags[[j]])
    }
    
    topic_frags_df = do.call(rbind.data.frame, df)
    topic_frags_df = topic_frags_df[,c("chr","start","end","RG")]
    topic_frags_df$counts  = rep(1,nrow(topic_frags_df))
    topic_frags_df$RG = paste0(colnames(topic_cells)[i],":",topic_frags_df$RG) 
    topic_frags_df$bc = do.call("rbind", strsplit(as.character(topic_frags_df$RG), "#"))[,2]
    topic_frags_df$lib = do.call("rbind", strsplit(as.character(topic_frags_df$RG), "#"))[,1]
    topic_frags_df$RG = paste0(topic_frags_df$lib,".",topic_frags_df$bc) 
    topic_frags_df$bc = NULL
    topic_frags_df$lib = NULL
    
    topic_frags_df <- topic_frags_df %>% arrange(chr, start, end)
    
    write.table(topic_frags_df,
                file = paste0("../topics/topic_", colnames(topic_cells)[i],
                             "_fragments.tsv"),
                sep="\t",
                quote=F,
                col.names = F, 
                row.names = F)
                
}

