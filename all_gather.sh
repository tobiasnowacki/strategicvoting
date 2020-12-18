#!/bin/sh
#SBATCH --time=1:00:00
#SBATCH --mem=20G
#SBATCH --partition=hns
#SBATCH --ntasks-per-node=6
#SBATCH --job-name="gather"
#SBATCH --error=gather.stderr
#SBATCH --output=gather.out
#SBATCH --mail-user=toby.nowacki@gmail.com
#SBATCH --mail-type=ALL

sbatch gather.sh 85 1
sbatch gather.sh 85 2
sbatch gather.sh 85 3
sbatch gather.sh 55 1
sbatch gather.sh 55 2
sbatch gather.sh 55 3
sbatch gather.sh 10 1
sbatch gather.sh 10 2
sbatch gather.sh 10 3
# Rscript code/new_iter_algo.R $1 $2