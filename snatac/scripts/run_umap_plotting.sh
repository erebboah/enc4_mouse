#!/bin/bash
#SBATCH --job-name=plot_archr   ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH -p standard               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1         ## number of cores the job needs
#SBATCH --output=plot_archr-%J.out ## output log file
#SBATCH --error=plot_archr-%J.err ## error log file
#SBATCH --mem=16G

source ~/miniconda3/bin/activate hpc3sc

Rscript ../scripts/plot_atac_umap.R
