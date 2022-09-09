#!/bin/bash
#SBATCH --job-name=archr_topics   ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH -p standard               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --cpus-per-task=8         ## number of cores the job needs
#SBATCH --output=archr_topics-%J.out ## output log file
#SBATCH --error=archr_topics-%J.err ## error log file
#SBATCH --mem=64G

source ~/miniconda3/bin/activate hpc3sc

mkdir ../topics
Rscript extract_topic_frags.R

cd ../topics

# zip and index fragment files
for f in *_fragments.tsv; do bgzip $f; done
for f in *_fragments.tsv.gz; do tabix -p bed $f; done
