library(tidyverse)

hats = read.delim("GO_term_summary_0004402.txt",row.names = NULL)
hats$Type = "HAT"
hats$GO_term_ID = "0004402"
hats$GO_term = "histone acetyltransferase activity"
hdacs = read.delim("GO_term_summary_0004407.txt",row.names = NULL)
hdacs$Type = "HDAC"
hdacs$GO_term_ID = "0004407"
hdacs$GO_term = "histone deacetylase activity"
hmts = read.delim("GO_term_summary_0042054.txt",row.names = NULL)
hmts$Type = "HMT"
hmts$GO_term_ID = "0042054"
hmts$GO_term = "histone methyltransferase activity"
hdms = read.delim("GO_term_summary_0032452.txt",row.names = NULL)
hdms$Type = "HDM"
hdms$GO_term_ID = "0032452"
hdms$GO_term = "histone demethylase activity"

med = read.delim("GO_term_summary_0016592.txt",row.names = NULL)
med$Type = "Mediator"
med$GO_term_ID = "0016592"
med$GO_term = "mediator complex"

tafs = read.delim("GO_term_summary_0006352.txt",row.names = NULL)
tafs$Type = "TAF"
tafs$GO_term_ID = "0006352"
tafs$GO_term = "DNA-templated transcription, initiation"


chrombind = read.delim("GO_term_summary_0003682.txt",row.names = NULL)
chrombind$Type = "Chromatin_binding"
chrombind$GO_term_ID = "0003682"
chrombind$GO_term = "chromatin binding"

chromorg = read.delim("GO_term_summary_0006325.txt",row.names = NULL)
chromorg$Type = "Chromatin_organization"
chromorg$GO_term_ID = "0006325"
chromorg$GO_term = "chromatin organization"

histone_regulators = rbind(unique(hats[,c("MGI.Gene.Marker.ID","Proteoform","Type","GO_term_ID","GO_term")]),
                          unique(hdacs[,c("MGI.Gene.Marker.ID","Proteoform","Type","GO_term_ID","GO_term")]),
                          unique(hmts[,c("MGI.Gene.Marker.ID","Proteoform","Type","GO_term_ID","GO_term")]),
                          unique(hdms[,c("MGI.Gene.Marker.ID","Proteoform","Type","GO_term_ID","GO_term")]),
                          unique(med[,c("MGI.Gene.Marker.ID","Proteoform","Type","GO_term_ID","GO_term")]),
                          unique(tafs[,c("MGI.Gene.Marker.ID","Proteoform","Type","GO_term_ID","GO_term")]),
                          unique(chrombind[,c("MGI.Gene.Marker.ID","Proteoform","Type","GO_term_ID","GO_term")]),
                          unique(chromorg[,c("MGI.Gene.Marker.ID","Proteoform","Type","GO_term_ID","GO_term")]))
colnames(histone_regulators) = c("gene_name","Evidence","Type","GO_term_ID","GO_term")

write.table(histone_regulators,file="histone_regulators.tsv",sep="\t",quote=F,row.names=F)

