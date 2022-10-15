#!/bin/bash
#SBATCH --job-name=topicproj    ## Name of the job
#SBATCH -A cosmos2022            ## account to charge 
#SBATCH -p standard               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --cpus-per-task=16         ## number of cores the job needs
#SBATCH --output=topproj-%J.out ## output log file
#SBATCH --error=topproj-%J.err ## error log file
#SBATCH --mem=128G

source ~/miniconda3/bin/activate hpc3sc

Rscript make_topic_archr_proj.R --tissue adrenal
Rscript make_topic_archr_proj.R --tissue gastrocnemius
Rscript make_topic_archr_proj.R --tissue heart
Rscript make_topic_archr_proj.R --tissue cortex
Rscript make_topic_archr_proj.R --tissue hippocampus
