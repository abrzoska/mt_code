#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=500G
#SBATCH --time=0-2
#SBATCH -p bigmem
#SBATCH -J maf_mapper
#SBATCH -o log.mapper
#SBATCH -A m2_jgu-nannospalax

module purge
/home/timlin/projects/mt_code/venv/bin/activate

python map_genes_to_cisregs.py