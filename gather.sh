#!/bin/sh
#SBATCH --time=3:00:00
#SBATCH --mem=200G
#SBATCH --partition=bigmem
#SBATCH --ntasks-per-node=6
#SBATCH --job-name="servercode"
#SBATCH --error=rep.stderr
#SBATCH --output=rep.out
#SBATCH --mail-user=toby.nowacki@gmail.com
#SBATCH --mail-type=ALL


echo $1
echo $2
module load R/4.0
module load cairo
Rscript code/gather_data.R $1 $2
# Rscript code/new_iter_algo.R $1 $2