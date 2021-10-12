#!/bin/sh
#SBATCH --time=2:00:00
#SBATCH --mem=64G
#SBATCH --partition=andyhall
#SBATCH --ntasks-per-node=24
#SBATCH --job-name="iter_four"
#SBATCH --error=scriptlog/rep4.stderr
#SBATCH --output=scriptlog/rep4.out
#SBATCH --mail-user=toby.nowacki@gmail.com
#SBATCH --mail-type=ALL

cd ..
echo $1
echo $2
module load R/4.0
Rscript code/new_iter_algo_four.R $1 $2
