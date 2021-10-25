#!/bin/env bash
#SBATCH -A single_cell_skeletal_muscle
#SBATCH -p fast
#SBATCH --mem 12G
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err

module load python/3.7
python analyse_sc.py --genelist Input/lama_etc_genes.tsv -o Output/MeanExpr_Lama_genes/
