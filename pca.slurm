#!/bin/bash
#SBATCH --job-name="pca"
#SBATCH --output="log/pca_%j.%N.out"
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=240Gb
#SBATCH --export=ALL
#SBATCH -t 48:00:00
#SBATCH --mail-user=adelorme@ucsd.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH -A csd403
source ~/.bashrc
cd ~/nemar/pca_averef
module load matlab
matlab -nodisplay -nosplash -nodesktop -r "compare_pca;"
