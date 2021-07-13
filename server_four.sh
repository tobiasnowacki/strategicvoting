#!/bin/sh
#SBATCH --time=4:00:00
#SBATCH --mem=266G
#SBATCH --partition=bigmem
#SBATCH --ntasks-per-node=4
#SBATCH --job-name="servercode"
#SBATCH --error=rep4.stderr
#SBATCH --output=rep4.out
#SBATCH --mail-user=toby.nowacki@gmail.com
#SBATCH --mail-type=ALL

echo $1
echo $2
module load R/4.0
Rscript code/new_iter_algo_four.R $1 $2
