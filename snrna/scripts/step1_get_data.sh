#!/bin/bash
#SBATCH --job-name=get_data    ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH -p standard               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --cpus-per-task=8         ## number of cores the job needs
#SBATCH --output=data-%J.out ## output log file
#SBATCH --error=data-%J.err ## error log file
#SBATCH --mem=128G

mkdir ../counts/
cd ../counts/

xargs -L 1 curl -O -J -L < ../ref/files.txt

for f in *.tar.gz
do
	mkdir "${f%.tar.gz}"
	tar xf "$f" -C "${f%.tar.gz}"
	rm "$f"
done

module load R
mkdir ../scrublet
mkdir ../seurat

Rscript ../scripts/merge_counts.R
