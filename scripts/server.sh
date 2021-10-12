#!/bin/sh
#SBATCH --time=2:00:00
#SBATCH --mem=64G
#SBATCH --partition=andyhall
#SBATCH --ntasks-per-node=24
#SBATCH --job-name="servercode"
#SBATCH --error=scriptlog/rep.stderr
#SBATCH --output=scriptlog/rep.out
#SBATCH --mail-user=toby.nowacki@gmail.com
#SBATCH --mail-type=ALL

cd ..
echo $1
echo $2
module load R/4.0
module load cairo
Rscript code/new_iter_algo.R $1 $2
