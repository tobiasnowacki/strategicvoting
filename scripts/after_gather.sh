#!/bin/sh
#SBATCH --time=1:00:00
#SBATCH --mem=20G
#SBATCH --partition=hns
#SBATCH --ntasks-per-node=6
#SBATCH --job-name="rest"
#SBATCH --error=scriptlog/rest
#SBATCH --output=scriptlog/rest
#SBATCH --mail-user=toby.nowacki@gmail.com
#SBATCH --mail-type=ALL

cd ..
Rscript code/new_convergence.R
Rscript code/new_summary.R
# Rscript code/new_iter_algo.R $1 $2