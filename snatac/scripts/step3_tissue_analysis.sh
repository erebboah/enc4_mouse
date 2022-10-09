#!/bin/bash
#SBATCH --job-name=tissues    ## Name of the job
#SBATCH -A cosmos2022            ## account to charge 
#SBATCH -p standard               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --cpus-per-task=16         ## number of cores the job needs
#SBATCH --output=archr_tissues-%J.out ## output log file
#SBATCH --error=archr_tissues-%J.err ## error log file
#SBATCH --mem=64G

source ~/miniconda3/bin/activate hpc3sc

Rscript archr_tissue_objects.R
