#!/bin/sh
#SBATCH --time=0:30:00
#SBATCH --mem=32G
#SBATCH --partition=hns
#SBATCH --ntasks-per-node=6
#SBATCH --job-name="servercode"
#SBATCH --error=rep_g.stderr
#SBATCH --output=rep_g.out
#SBATCH --mail-user=toby.nowacki@gmail.com
#SBATCH --mail-type=ALL


echo $1
echo $2
module load R/4.0
module load cairo
Rscript code/gather_summary_four.R $1 $2
# Rscript code/new_iter_algo.R $1 $2
