suppressPackageStartupMessages(library(ArchR))
library(parallel)
suppressPackageStartupMessages(library(Seurat))
library(reshape2)
library(Repitools)
library(tidyverse)
addArchRThreads(threads = 16)
addArchRGenome("mm10")
library(optparse)
set.seed(1234)

option_list = list(
  make_option(c("-t", "--tissue"), action="store", default=NA, type='character',
              help="hippocampus, cortex, heart, adrenal, or gastrocnemius"))

opt = parse_args(OptionParser(option_list=option_list))
tissue = opt$tissue

setwd("../archr/")

projectname = paste0("ENC4_Mouse_",str_to_title(tissue))
proj = loadArchRProject(path = projectname)
meta = getCellColData(proj)

fn = paste0("../ref/",tolower(tissue),"_selectedCells_top0.05_0.05min_score_50min_cells.csv")
topic_cells = read.csv(fn,row.names=1)

if (tolower(tissue) == "heart"){
    abrv = "HT"
    } else if (tolower(tissue) == "gastrocnemius"){
    abrv = "GC"
    } else if (tolower(tissue) == "hippocampus"){
    abrv = "HC"   
    } else if (tolower(tissue) == "cortex"){
    abrv = "CX" 
    } else if (tolower(tissue) == "adrenal"){
    abrv = "AD" 
}

colnames(topic_cells) = gsub("Topic_",abrv,colnames(topic_cells))

for (i in 1:ncol(topic_cells)) {
    this_topic = topic_cells[,i]
    this_topic = this_topic[this_topic != ""] # no blanks
    this_topic = do.call("rbind", strsplit(as.character(this_topic), "-1"))[,1]
    print(table(this_topic %in% meta$cellID))
    
    meta_subset = meta[meta$cellID %in% this_topic,]
    meta_subset = meta_subset[match(this_topic, meta_subset$cellID),]
    
    topic_frags = getFragmentsFromProject(ArchRProj = proj, cellNames = rownames(meta_subset))
  
    df = list()
    for (j in 1:length(topic_frags)){
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
    
    fn = paste0("../archr/",projectname,"_Topics/",colnames(topic_cells)[i],
                             "_fragments.tsv")
    print(fn)
    
    write.table(topic_frags_df,
                file = paste0("../archr/",projectname,"_Topics/",colnames(topic_cells)[i],
                             "_fragments.tsv"),
                sep="\t",
                quote=F,
                col.names = F, 
                row.names = F)
    
}