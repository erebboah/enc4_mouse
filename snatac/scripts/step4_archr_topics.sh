#!/bin/bash
#SBATCH --job-name=topics_archr   ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH -p standard               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --cpus-per-task=16         ## number of cores the job needs
#SBATCH --output=topics_archr-%J.out ## output log file
#SBATCH --error=topics_archr-%J.err ## error log file
#SBATCH --mem=64G

source ~/miniconda3/bin/activate hpc3sc

Rscript make_topic_archr_proj.R
