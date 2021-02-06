#!/bin/sh
#SBATCH --time=1:00:00
#SBATCH --mem=20G
#SBATCH --partition=hns
#SBATCH --ntasks-per-node=6
#SBATCH --job-name="randproc"
#SBATCH --error=randproc.stderr
#SBATCH --output=randproc.out
#SBATCH --mail-user=toby.nowacki@gmail.com
#SBATCH --mail-type=ALL

Rscript code/new_random_process.R
# Rscript code/new_iter_algo.R $1 $2