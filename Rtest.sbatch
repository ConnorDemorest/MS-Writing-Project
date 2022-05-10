#!/bin/bash
#SBATCH --partition classpart
#SBATCH --account classroom
#SBATCH --cpus-per-task 40
#SBATCH --mem 40000
#SBATCH --time 0-00:01:00
#SBATCH --job-name Rtest
#SBATCH --output Rtest-%j.out
#SBATCH --error Rtest-%j.err
#SBATCH --mail-type ALL
module load Anaconda3/3.5.1
source activate $HOME/renv/
Rscript MHT_HPC_script.R
