suppressPackageStartupMessages(library(ArchR))
library(parallel)
library(viridis)
suppressPackageStartupMessages(library(Seurat))
set.seed(1234)
library(reshape2)
addArchRThreads(threads = 1) 
addArchRGenome("mm10")

setwd("/share/crsp/lab/seyedam/share/enc4_mouse/snatac/archr/")

proj = loadArchRProject(path = "ENC4_Mouse/")


topic_cells = read.csv("/share/crsp/lab/seyedam/share/enc4_mouse/snatac/ref/selectedcells_top5%_50cells.csv",row.names=1)

pathToMacs2 <- findMacs2()
meta = getCellColData(proj)



combined.peaks = readRDS("/share/crsp/lab/seyedam/share/enc4_mouse/snatac/archr/topics/all_combined_peaks.rds")
                         
cell_scores = list()
marker_list = list()

for (i in 1:ncol(topic_cells)){
    meta_subset = meta[meta$rna_bc %in% unique(topic_cells[,i][topic_cells[,i] != ""]),]
    archr_topic = subsetCells(ArchRProj = proj, cellNames = rownames(meta_subset))
    archr_topic$Topic = rep(colnames(topic_cells)[i],nrow(meta_subset))
    
    archr_topic = addPeakSet(ArchRProj = archr_topic,
                         peakSet = combined.peaks,
                         force = TRUE)
    
    archr_topic <- addPeakMatrix(archr_topic) # quantify
    
    
    rangedsumexpt = getMatrixFromProject(ArchRProj = archr_topic,
                                            useMatrix = "PeakMatrix",
                                            verbose = TRUE,
                                            binarize = FALSE) # binarize?
    cell_scores[[i]] = assay(rangedsumexpt)

}

cell_score_avg <- data.frame(matrix(ncol = ncol(topic_cells), nrow = length(combined.peaks)))
colnames(cell_score_avg) <- colnames(topic_cells)

for (i in 1:ncol(topic_cells)) {
    cell_score_avg[,i] = rowMeans(cell_scores[[i]])
}

saveRDS(cell_score_avg,"/share/crsp/lab/seyedam/share/enc4_mouse/snatac/archr/topics/topics_cell_score_matrix.rds")


