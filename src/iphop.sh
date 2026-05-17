#!/bin/bash
#SBATCH --job-name=iphop
#SBATCH --time=48:00:00
#SBATCH --partition=small
#SBATCH --account=project_2001499
#SBATCH --mem=50G
#SBATCH --cpus-per-task=12
#SBATCH --gres=nvme:100

export PATH=/projappl/project_2001499/iphop/bin:$PATH

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

iphop predict --fa_file /scratch/project_2001499/$USER/MBDP_Metagenomics_2026/05_VIROMICS/vOTUs/vOTUs.fna \
--min_score 75 \
--db_dir /scratch/project_2001499/DBs/IPHOP_Jun_2025_pub_rw \
--out_dir /scratch/project_2001499/$USER/MBDP_Metagenomics_2026/05_VIROMICS/IPHOP \
-t $SLURM_CPUS_PER_TASK \
--single_thread_wish
