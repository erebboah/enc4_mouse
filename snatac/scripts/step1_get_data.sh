#!/bin/bash
#SBATCH --job-name=get_data  ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH -p standard               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --cpus-per-task=8         ## number of cores the job needs
#SBATCH --output=atacdata-%J.out ## output log file
#SBATCH --error=atacdata-%J.err ## error log file
#SBATCH --mem=128G

mkdir ../fragments/
cd ../fragments/

xargs -L 1 curl -O -J -L < ../ref/encode_10x_fragment_files.txt

for f in *.tar.gz
do
	mkdir "${f%.tar.gz}"
	tar xf "$f" -C "${f%.tar.gz}"
	rm "$f"
done

cd ../scripts/
source ~/miniconda3/bin/activate hpc3sc

mkdir ../archr
Rscript make_archr_proj.R
