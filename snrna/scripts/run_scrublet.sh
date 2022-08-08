#!/bin/bash
#SBATCH --job-name=pscrub    ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH -p standard               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --cpus-per-task=8         ## number of cores the job needs
#SBATCH --output=pscrub-%J.out ## output log file
#SBATCH --error=pscrub-%J.err ## error log file
#SBATCH --mem=64G
#SBATCH --time 3-00:00

source /data/homezvol2/erebboah/miniconda3/bin/activate hpc3sc

mkdir ../scrublet
mkdir ../seurat

python3 run_scrublet.py
