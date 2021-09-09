#!/bin/sh
<<<<<<< HEAD
<<<<<<< HEAD
#SBATCH --time=3:00:00
#SBATCH --mem=200G
#SBATCH --partition=bigmem
#SBATCH --ntasks-per-node=6
#SBATCH --job-name="servercode"
#SBATCH --error=rep.stderr
#SBATCH --output=rep.out
=======
=======
>>>>>>> 95425d540f4efbe579af50da3c43f21f50cc9383
#SBATCH --time=2:00:00
#SBATCH --mem=64G
#SBATCH --partition=andyhall
#SBATCH --ntasks-per-node=6
#SBATCH --job-name="gather"
#SBATCH --error=g.stderr
#SBATCH --output=g.out
<<<<<<< HEAD
>>>>>>> 95425d540f4efbe579af50da3c43f21f50cc9383
=======
>>>>>>> 95425d540f4efbe579af50da3c43f21f50cc9383
#SBATCH --mail-user=toby.nowacki@gmail.com
#SBATCH --mail-type=ALL


echo $1
echo $2
module load R/4.0
module load cairo
Rscript code/new_gather_data.R $1 $2
# Rscript code/new_iter_algo.R $1 $2
