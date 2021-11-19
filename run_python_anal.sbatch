#!/bin/env bash
#SBATCH -A single_cell_skeletal_muscle
#SBATCH -p fast
#SBATCH --cpus-per-task 20
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err

module load python/3.7

data_dir="/projects/single_cell_skeletal_muscle/"
genelist=${data_dir}/Input/my_genes_of_interest.tsv
metadata=${data_dir}/GSE_repo/Metadata.txt
datafile=${data_dir}/GSE_repo/normalised_expression_matrix.txt
outdir=${data_dir}/Output/

module load python/3.7

python analyse_sc.py -g ${genelist} -d ${datafile} -m ${metadata} -o ${outdir}
