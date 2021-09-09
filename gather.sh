#!/bin/sh
#SBATCH --time=1:00:00
#SBATCH --mem=64G
#SBATCH --partition=andyhall
#SBATCH --ntasks-per-node=6
#SBATCH --job-name="gather"
#SBATCH --error=g.stderr
#SBATCH --output=g.out
#SBATCH --mail-user=toby.nowacki@gmail.com
#SBATCH --mail-type=ALL

echo $1
echo $2
module load R/4.0
module load cairo
Rscript code/new_gather_data.R $1 $2
# Rscript code/new_iter_algo.R $1 $2
