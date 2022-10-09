#!/bin/bash
#SBATCH --job-name=submit    ## Name of the job
#SBATCH -A cosmos2022           ## account to charge 
#SBATCH -p standard               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --cpus-per-task=8         ## number of cores the job needs
#SBATCH --output=files-%J.out ## output log file
#SBATCH --error=files-%J.err ## error log file
#SBATCH --mem=256G
#SBATCH --time 3-00:00

source ~/miniconda3/bin/activate hpc3sc

Rscript Submission_files.R
