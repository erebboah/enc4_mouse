#!/bin/bash
#SBATCH --job-name=tissue_int   ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH -p highmem               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --cpus-per-task=8         ## number of cores the job needs
#SBATCH --output=int_tissue-%J.out ## output log file
#SBATCH --error=int_tissue-%J.err ## error log file
#SBATCH --mem=350G

module load R

Rscript integrate_gene_subset.R
