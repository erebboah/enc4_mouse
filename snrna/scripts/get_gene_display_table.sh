#!/bin/bash
#SBATCH --job-name=get_table    ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH -p standard               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --cpus-per-task=16         ## number of cores the job needs
#SBATCH --output=genes-%J.out ## output log file
#SBATCH --error=genes-%J.err ## error log file
#SBATCH --mem=16G

source /data/homezvol2/erebboah/miniconda3/bin/activate hpc3sc

cd ../ref

wget "http://www.informatics.jax.org/go/report.txt?goID=GO:0004402&results=83&startIndex=0&sort=term&dir="
mv report.txt* GO_term_summary_0004402.txt

wget "http://www.informatics.jax.org/go/report.txt?goID=GO:0004407&results=70&startIndex=0&sort=term&dir="
mv report.txt* GO_term_summary_0004407.txt

wget "http://www.informatics.jax.org/go/report.txt?goID=GO:0042054&results=173&startIndex=0&sort=term&dir="
mv report.txt* GO_term_summary_0042054.txt

wget "http://www.informatics.jax.org/go/report.txt?goID=GO:0032452&results=113&startIndex=0&sort=term&dir="
mv report.txt* GO_term_summary_0032452.txt

wget "http://www.informatics.jax.org/go/report.txt?goID=GO:0016592&results=97&startIndex=0&sort=term&dir="
mv report.txt* GO_term_summary_0016592.txt

wget "http://www.informatics.jax.org/go/report.txt?goID=GO:0006352&results=247&startIndex=0&sort=term&dir="
mv report.txt* GO_term_summary_0006352.txt

wget "http://www.informatics.jax.org/go/report.txt?goID=GO:0003682&results=1209&startIndex=0&sort=term&dir="
mv report.txt* GO_term_summary_0003682.txt

wget "http://www.informatics.jax.org/go/report.txt?goID=GO:0006325&results=1059&startIndex=0&sort=term&dir="
mv report.txt* GO_term_summary_0006325.txt

Rscript ../scripts/get_gene_display_table.R
