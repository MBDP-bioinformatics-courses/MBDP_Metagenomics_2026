#!/bin/bash
#SBATCH --job-name=gtdbtk
#SBATCH --account=project_2001499
#SBATCH --time=6:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --gres=nvme:48
#SBATCH --output=%x-%j.out

# initialise
module load gtdbtk

# you need to add below the path to the folder containing all the MAGs 
# e.g. /scratch/project_2001499/$USER/MBDP_Metagenomics_2026/06_ANVIO/ALL_BINS
genome_dir=

# run gtdbtk
gtdbtk \
  classify_wf \
  --genome_dir $genome_dir \
  --out_dir GTDB \
  --extension fa \
  --cpus $SLURM_CPUS_PER_TASK \
  --pplacer_cpus 2
