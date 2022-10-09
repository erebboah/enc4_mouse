#!/bin/bash
#SBATCH --job-name=archrproj  ## Name of the job
#SBATCH -A cosmos2022            ## account to charge 
#SBATCH -p standard               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --cpus-per-task=8         ## number of cores the job needs
#SBATCH --output=archrproj-%J.out ## output log file
#SBATCH --error=archrproj-%J.err ## error log file
#SBATCH --mem=128G

source ~/miniconda3/bin/activate hpc3sc

Rscript process_archr_proj.R
