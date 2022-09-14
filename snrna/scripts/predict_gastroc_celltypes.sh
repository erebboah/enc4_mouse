#!/bin/bash
#SBATCH --job-name=gc_pred    ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH -p standard               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --cpus-per-task=8         ## number of cores the job needs
#SBATCH --output=gc-%J.out ## output log file
#SBATCH --error=gc-%J.err ## error log file
#SBATCH --mem=128G

source /data/homezvol2/erebboah/miniconda3/bin/activate hpc3sc

mkdir ../plots/gastrocnemius
mkdir ../plots/gastrocnemius/annotation
Rscript predict_gastroc_celltypes.R
