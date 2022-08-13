#!/bin/bash
#SBATCH --job-name=arraycellbend    ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH -p standard               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --cpus-per-task=8         ## number of cores the job needs
#SBATCH --array=1-40               ## number of tasks to launch (number of samples)
#SBATCH --output=arraycellbend-%J.out ## output log file
#SBATCH --error=arraycellbend-%J.err ## error log file
#SBATCH --mem=128G
#SBATCH --time 14-00:00

source /data/homezvol2/erebboah/miniconda3/bin/activate cellbender2

file="../ref/10x_cellbender_settings.csv"

settings=`cat $file | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`

infile=$(echo $settings | cut -d ',' -f 1)
outfile=$infile"/cellbender.h5"
cells=$(echo $settings | cut -d ',' -f 2)
total_drops=$(echo $settings | cut -d ',' -f 3)

# run cellbender
cellbender remove-background \
	--input $infile \
	--output $outfile \
	--expected-cells $cells \
	--total-droplets-included $total_drops \
	--epochs 150 \
	--fpr 0.01
