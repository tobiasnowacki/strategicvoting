#!/bin/sh
#SBATCH --time=2:00:00
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
#SBATCH --mem=194G
#SBATCH --partition=bigmem
#SBATCH --ntasks-per-node=4
=======
#SBATCH --mem=64G
#SBATCH --partition=andyhall
#SBATCH --ntasks-per-node=24
>>>>>>> 95425d540f4efbe579af50da3c43f21f50cc9383
=======
#SBATCH --mem=64G
#SBATCH --partition=andyhall
#SBATCH --ntasks-per-node=24
>>>>>>> 95425d540f4efbe579af50da3c43f21f50cc9383
=======
#SBATCH --mem=64G
#SBATCH --partition=andyhall
#SBATCH --ntasks-per-node=24
>>>>>>> 95425d540f4efbe579af50da3c43f21f50cc9383
#SBATCH --job-name="servercode"
#SBATCH --error=rep.stderr
#SBATCH --output=rep.out
#SBATCH --mail-user=toby.nowacki@gmail.com
#SBATCH --mail-type=ALL

echo $1
echo $2
module load R/4.0
module load cairo
Rscript code/new_iter_algo.R $1 $2
