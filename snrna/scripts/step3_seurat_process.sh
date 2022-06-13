#!/bin/bash
#SBATCH --job-name=seurat_process   ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH -p standard               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --cpus-per-task=8         ## number of cores the job needs
#SBATCH --output=filt-%J.out ## output log file
#SBATCH --error=filt-%J.err ## error log file
#SBATCH --mem=256G

module load R

Rscript seurat_process.R
