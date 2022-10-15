#!/bin/bash
#SBATCH --job-name=topicfrags    ## Name of the job
#SBATCH -A cosmos2022            ## account to charge 
#SBATCH -p standard               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --cpus-per-task=16         ## number of cores the job needs
#SBATCH --output=frags-%J.out ## output log file
#SBATCH --error=frags-%J.err ## error log file
#SBATCH --mem=128G

source ~/miniconda3/bin/activate hpc3sc

wget -O ../ref/adrenal_selectedCells_top0.05_0.05min_score_50min_cells.csv https://raw.githubusercontent.com/nargesr/Topics-modeling/main/TFs_mirhgs_chromreg_reproducible_final/Adrenal_parse_10x_20/selectedCells_top0.05_0.05min_score_50min_cells.csv?token=GHSAT0AAAAAABYO75XLE4AL5M2MHQHFUT2OY2ITAQQ

wget -O ../ref/cortex_selectedCells_top0.05_0.05min_score_50min_cells.csv https://raw.githubusercontent.com/nargesr/Topics-modeling/main/TFs_mirhgs_chromreg_reproducible_final/Cortex_parse_10x_20/selectedCells_top0.05_0.05min_score_50min_cells.csv?token=GHSAT0AAAAAABYO75XK3PPFRHZEAZCZFXWOY2ITB3Q

wget -O ../ref/gastrocnemius_selectedCells_top0.05_0.05min_score_50min_cells.csv https://raw.githubusercontent.com/nargesr/Topics-modeling/main/TFs_mirhgs_chromreg_reproducible_final/Gastrocnemius_parse_10x_25/selectedCells_top0.05_0.05min_score_50min_cells.csv?token=GHSAT0AAAAAABYO75XLDPHRIJ7AZ2IQF4YIY2ITCEA

wget -O ../ref/heart_selectedCells_top0.05_0.05min_score_50min_cells.csv https://raw.githubusercontent.com/nargesr/Topics-modeling/main/TFs_mirhgs_chromreg_reproducible_final/Heart_parse_10x_20/selectedCells_top0.05_0.05min_score_50min_cells.csv?token=GHSAT0AAAAAABYO75XL3TN3BRKIQAZE2INUY2ITCMA

wget -O ../ref/hippocampus_selectedCells_top0.05_0.05min_score_50min_cells.csv https://raw.githubusercontent.com/nargesr/Topics-modeling/main/TFs_mirhgs_chromreg_reproducible_final/Hippocampus_parse_10x_25/selectedCells_top0.05_0.05min_score_50min_cells.csv?token=GHSAT0AAAAAABYO75XKFZPCJ4JTG272VCLYY2ITCUA

mkdir ../archr/ENC4_Mouse_Adrenal_Topics
mkdir ../archr/ENC4_Mouse_Cortex_Topics
mkdir ../archr/ENC4_Mouse_Gastrocnemius_Topics
mkdir ../archr/ENC4_Mouse_Hippocampus_Topics
mkdir ../archr/ENC4_Mouse_Heart_Topics

Rscript extract_topic_frags.R --tissue adrenal
cd ../archr/ENC4_Mouse_Adrenal_Topics
for f in *_fragments.tsv; do bgzip $f; done
for f in *_fragments.tsv.gz; do tabix -p bed $f; done

cd ../../scripts
Rscript extract_topic_frags.R --tissue cortex 
cd ../archr/ENC4_Mouse_Cortex_Topics
for f in *_fragments.tsv; do bgzip $f; done
for f in *_fragments.tsv.gz; do tabix -p bed $f; done

cd ../../scripts
Rscript extract_topic_frags.R --tissue heart
cd ../archr/ENC4_Mouse_Heart_Topics
for f in *_fragments.tsv; do bgzip $f; done
for f in *_fragments.tsv.gz; do tabix -p bed $f; done

cd ../../scripts
Rscript extract_topic_frags.R --tissue hippocampus
cd ../archr/ENC4_Mouse_Hippocampus_Topics
for f in *_fragments.tsv; do bgzip $f; done
for f in *_fragments.tsv.gz; do tabix -p bed $f; done

cd ../../scripts
Rscript extract_topic_frags.R --tissue gastrocnemius
cd ../archr/ENC4_Mouse_Gastrocnemius_Topics
for f in *_fragments.tsv; do bgzip $f; done
for f in *_fragments.tsv.gz; do tabix -p bed $f; done
