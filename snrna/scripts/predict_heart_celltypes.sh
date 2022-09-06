#!/bin/bash
#SBATCH --job-name=ht_pred    ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH -p standard               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --cpus-per-task=16         ## number of cores the job needs
#SBATCH --output=ht2-%J.out ## output log file
#SBATCH --error=ht2-%J.err ## error log file
#SBATCH --mem=256G

source /data/homezvol2/erebboah/miniconda3/bin/activate hpc3sc

cd ../ref/external_data
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8810/E-MTAB-8810.processed.1.zip

unzip E-MTAB-8810.processed.1.zip

wget https://cellgeni.cog.sanger.ac.uk/heartcellatlas/data/global_raw.h5ad

cd ../../scripts

mkdir ../plots/heart
mkdir ../plots/heart/annotation
Rscript predict_heart_celltypes.R
