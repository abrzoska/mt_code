#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=100MB
#SBATCH -t 15
#SBATCH -p smp
#SBATCH -J maf_mapper
#SBATCH -o mapper_script.log
#SBATCH -A m2_jgu-nannospalax

module purge
module load lang/SciPy-bundle

source ~/venv/bin/activate

python /home/timlin/projects/mt_code/test_script_memory_profiler.py