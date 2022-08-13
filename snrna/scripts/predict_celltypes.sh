#!/bin/bash
#SBATCH --job-name=hc_prd    ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH -p standard               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --cpus-per-task=8         ## number of cores the job needs
#SBATCH --output=hc2-%J.out ## output log file
#SBATCH --error=hc2-%J.err ## error log file
#SBATCH --mem=128G

source /data/homezvol2/erebboah/miniconda3/bin/activate hpc3sc

Rscript predict_celltypes.R
