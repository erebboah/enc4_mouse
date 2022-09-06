#!/bin/bash
#SBATCH --job-name=scrub    ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH -p standard               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --cpus-per-task=8         ## number of cores the job needs
#SBATCH --output=scrub-%J.out ## output log file
#SBATCH --error=scrub-%J.err ## error log file
#SBATCH --mem=128G
#SBATCH --time 3-00:00

source ~/miniconda3/bin/activate hpc3sc

#Rscript format_cellbender_output.R
python3 run_scrublet.py

mkdir ../plots
