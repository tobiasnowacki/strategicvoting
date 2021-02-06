#!/bin/sh
#SBATCH --time=1:00:00
#SBATCH --mem=200G
#SBATCH --partition=bigmem
#SBATCH --ntasks-per-node=6
#SBATCH --job-name="conj"
#SBATCH --error=conj.stderr
#SBATCH --output=conj.out
#SBATCH --mail-user=toby.nowacki@gmail.com
#SBATCH --mail-type=ALL

module load R/4.0
module load cairo
Rscript code/new_conjecture.R 
# Rscript code/new_iter_algo.R $1 $2
