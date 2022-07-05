#!/bin/bash
#SBATCH -c 32
#SBATCH --mem=96GB
#SBATCH -t 24:00:00
#SBATCH -p parallel
#SBATCH -J maf_mapper_indel3
#SBATCH -o mapper_script.log
#SBATCH -A m2_jgu-nannospalax

module purge
module load lang/SciPy-bundle

source ~/venv/bin/activate

python /home/timlin/projects/mt_code/run_streamlined_genome_analysis.py 3
