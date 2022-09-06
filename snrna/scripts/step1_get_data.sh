#!/bin/bash
#SBATCH --job-name=get_data    ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH -p standard               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --cpus-per-task=16         ## number of cores the job needs
#SBATCH --output=data-%J.out ## output log file
#SBATCH --error=data-%J.err ## error log file
#SBATCH --mem=256G

source ~/miniconda3/bin/activate hpc3sc

mkdir ../counts_parse/
cd ../counts_parse/

xargs -L 1 curl -O -J -L < ../ref/encode_parse_scrna_files.txt

for f in *.tar.gz
do
	mkdir "${f%.tar.gz}"
	tar xf "$f" -C "${f%.tar.gz}"
	rm "$f"
done

rm metadata.tsv
cd ../scripts

mkdir ../counts_10x/
cd ../counts_10x/

xargs -L 1 curl -O -J -L < ../ref/encode_10x_scrna_files.txt

for f in *.tar.gz
do
	mkdir "${f%.tar.gz}"
	tar xf "$f" -C "${f%.tar.gz}"
	rm "$f"
done

rm metadata.tsv

cd ../scripts

mkdir ../scrublet/
Rscript get_counts_parse.R
Rscript get_counts_10x.R
