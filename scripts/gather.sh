#!/bin/sh
#SBATCH --time=0:30:00
#SBATCH --mem=32G
#SBATCH --partition=andyhall
#SBATCH --ntasks-per-node=1
#SBATCH --job-name="gather"
#SBATCH --error=scriptlog/g.stderr
#SBATCH --output=scriptlog/g.out
#SBATCH --mail-user=toby.nowacki@gmail.com
#SBATCH --mail-type=ALL

echo $1
echo $2
module load R/4.0
module load cairo
Rscript ../code/new_gather_data.R $1 $2
# Rscript code/new_iter_algo.R $1 $2
