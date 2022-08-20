#!/bin/bash
#SBATCH --job-name=gc_int    ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH -p standard               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --cpus-per-task=8         ## number of cores the job needs
#SBATCH --output=gc1-%J.out ## output log file
#SBATCH --error=gc1-%J.err ## error log file
#SBATCH --mem=128G

source /data/homezvol2/erebboah/miniconda3/bin/activate hpc3sc

mkdir ../seurat
Rscript integrate_parse_10x.R --tissue gastrocnemius
