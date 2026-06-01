#!/bin/bash
#SBATCH --job-name=anvio
#SBATCH --account=project_2001499
#SBATCH --time=6:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=20
#SBATCH --gres=nvme:48
#SBATCH --output=%x-%j.out

# initialise
export PATH="/projappl/project_2001499/anvio-dev/bin/:/projappl/project_2001499/anvio-gh/:$PATH"
export PYTHONPATH="/projappl/project_2001499/anvio-gh/:$PYTHONPATH"

# first assembly
anvi-run-kegg-kofams \
  -c path-to-CONTIGS.db \
  -T $SLURM_CPUS_PER_TASK \
  --kegg-data-dir /scratch/project_2001499/DBs/anvio_kegg_data

# second assembly
anvi-run-kegg-kofams \
  -c path-to-CONTIGS.db \
  -T $SLURM_CPUS_PER_TASK \
  --kegg-data-dir /scratch/project_2001499/DBs/anvio_kegg_data
