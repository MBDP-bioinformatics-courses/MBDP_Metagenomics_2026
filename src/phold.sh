#!/bin/bash
#SBATCH --job-name=phold
#SBATCH --account=project_2001499
#SBATCH --time=04:00:00
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=10G
#SBATCH --gres=gpu:v100:1

export PATH="/projappl/project_2001499/phold/bin:$PATH"

cd /scratch/project_2001499/$USER/MBDP_Metagenomics_2026/05_VIROMICS

phold run \
-i /scratch/project_2001499/$USER/MBDP_Metagenomics_2026/05_VIROMICS/vOTUs/vOTUs.fna \
-o PHOLD \
-d /scratch/project_2001499/DBs/PHOLD_DB \
-f &> /scratch/project_2001499/$USER/MBDP_Metagenomics_2026/00_LOGS/phold.log
