#!/bin/sh
#SBATCH --time=1:00:00
#SBATCH --mem=20G
#SBATCH --partition=hns
#SBATCH --ntasks-per-node=6
#SBATCH --job-name="gather"
#SBATCH --error=scriptlog/gather.stderr
#SBATCH --output=scriptlog/gather.out
#SBATCH --mail-user=toby.nowacki@gmail.com
#SBATCH --mail-type=ALL

module load R/4.0
module load cairo
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
