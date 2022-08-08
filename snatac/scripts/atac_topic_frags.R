suppressPackageStartupMessages(library(ArchR))
library(parallel)
library(viridis)
suppressPackageStartupMessages(library(Seurat))
set.seed(1234)
library(reshape2)
library(Repitools)
addArchRThreads(threads = 1) 
addArchRGenome("mm10")


setwd("/share/crsp/lab/seyedam/share/enc4_mouse/snatac/archr/")

proj = loadArchRProject(path = "ENC4_Mouse/")


topic_cells = read.csv("/share/crsp/lab/seyedam/share/enc4_mouse/snatac/ref/selectedcells_top5%_50cells.csv",row.names=1)
meta = getCellColData(proj)


topic_frags = list()
topic_frags_df = list()
for (i in 75:ncol(topic_cells)) {
    print(table(unique(topic_cells[,i][topic_cells[,i] != ""]) %in% meta$rna_bc))
    
    meta_subset = meta[meta$rna_bc %in% unique(topic_cells[,i][topic_cells[,i] != ""]),]
            
    topic_frags[[i]] = getFragmentsFromProject(ArchRProj = proj,
                            cellNames = rownames(meta_subset))
    df = list()
    for (j in 1:40){
        df[[j]] <- annoGR2DF(topic_frags[[i]][[j]])
    }
    
    topic_frags_df[[i]] = do.call(rbind.data.frame, df)
    topic_frags_df[[i]] = topic_frags_df[[i]][,c("chr","start","end","RG")]
    topic_frags_df[[i]]$counts  = rep(1,nrow(topic_frags_df[[i]]))
    
    write.table(topic_frags_df[[i]],
                file = paste0("../topics/topic_", colnames(topic_cells)[i],
                             "_fragments.tsv"),
                sep="\t",
                quote=F,
                col.names = F, 
                row.names = F)
}


