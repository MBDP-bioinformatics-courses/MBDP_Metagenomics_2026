#!/bin/bash
#SBATCH --job-name=gm
#SBATCH --time=06:00:00
#SBATCH --partition=small
#SBATCH --account=project_2001499
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4
#SBATCH --gres=nvme:50

export PATH="/projappl/project_2001499/genomad/bin:$PATH" 

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

while read i
do
genomad end-to-end \
--cleanup \
--splits 16 \
/scratch/project_2001499/$USER/MBDP_Metagenomics_2026/03_ASSEMBLY/${i}_flye/assembly.fasta \
/scratch/project_2001499/$USER/MBDP_Metagenomics_2026/05_VIROMICS/GENOMAD/${i} \
/scratch/project_2001499/DBs/genomad_db \
--threads $SLURM_CPUS_PER_TASK &> /scratch/project_2001499/$USER/MBDP_Metagenomics_2026/00_LOGS/genomad_${i}.log
done < $1
