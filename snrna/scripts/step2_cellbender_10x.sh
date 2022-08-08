#!/bin/bash
#SBATCH --job-name=cellbend    ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH -p standard               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --cpus-per-task=8         ## number of cores the job needs
#SBATCH --output=cellbend-%J.out ## output log file
#SBATCH --error=cellbend-%J.err ## error log file
#SBATCH --mem=128G
#SBATCH --time 14-00:00

source /data/homezvol2/erebboah/miniconda3/bin/activate cellbender2

for d in /share/crsp/lab/seyedam/share/enc4_mouse/snrna/counts_10x/*/;
do

infile=$d
outfile=$d"cellbender.h5"
settings=$(head -n 1 $infile"cellbender_options.csv")
cells=$(echo $settings | cut -d ',' -f 2)
total_drops=$(echo $settings | cut -d ',' -f 1)

# run cellbender
cellbender remove-background \
	--input $infile \
	--output $outfile \
	--expected-cells $cells \
	--total-droplets-included $total_drops \
	--epochs 150 \
	--fpr 0.01
done
