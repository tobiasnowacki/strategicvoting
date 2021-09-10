#!/bin/sh
#SBATCH --time=2:00:00
#SBATCH --mem=200G
#SBATCH --partition=bigmem
#SBATCH --ntasks-per-node=6
#SBATCH --job-name="randomiter"
#SBATCH --error=scriptlog/randomiter.stderr
#SBATCH --output=scriptlog/randomiter.out
#SBATCH --mail-user=toby.nowacki@gmail.com
#SBATCH --mail-type=ALL


echo $1
module load R/4.0
module load cairo
Rscript ../code/new_random_start.R $1
# Rscript code/new_iter_algo.R $1 $2
