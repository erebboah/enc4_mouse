#!/bin/bash
#SBATCH --job-name=brainref    ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH -p standard               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --cpus-per-task=16         ## number of cores the job needs
#SBATCH --output=brainref-%J.out ## output log file
#SBATCH --error=brainref-%J.err ## error log file
#SBATCH --mem=256G

mkdir ref/external_data
cd ref/external_data
wget https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hpf_10x/matrix.csv
wget https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hpf_10x/metadata.csv

sed -n -e '1,100000p' matrix.csv > matrix_1.csv
sed -n -e '100001,200000p' matrix.csv > matrix_2.csv
sed -n -e '200001,300000p' matrix.csv > matrix_3.csv
sed -n -e '300001,400000p' matrix.csv > matrix_4.csv
sed -n -e '400001,500000p' matrix.csv > matrix_5.csv
sed -n -e '500001,600000p' matrix.csv > matrix_6.csv
sed -n -e '600001,700000p' matrix.csv > matrix_7.csv
sed -n -e '700001,800000p' matrix.csv > matrix_8.csv
sed -n -e '800001,900000p' matrix.csv > matrix_9.csv
sed -n -e '900001,1000000p' matrix.csv > matrix_10.csv
sed -n -e '1000001,1100000p' matrix.csv > matrix_11.csv
sed -n -e '1100001,$p' matrix.csv > matrix_12.csv

rm matrix.csv

Rscript ../../subsample_brain_data.R
