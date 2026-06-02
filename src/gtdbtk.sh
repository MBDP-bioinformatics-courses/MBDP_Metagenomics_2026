#!/bin/bash
#SBATCH --job-name=gtdbtk
#SBATCH --account=project_2001499
#SBATCH --time=4:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=16
#SBATCH --gres=nvme:48
#SBATCH --output=%x-%j.out

# initialise
module load gtdbtk

# you need to add these below:
# summary_folder_1: path to the summary folder of assembly #1 (e.g. /scratch/project_2001499/$USER/MBDP_Metagenomics_2026/06_ANVIO/ERR5000342/METABAT_SUMMARY)
# summary_folder_2: path to the summary folder of assembly #2 (e.g. /scratch/project_2001499/$USER/MBDP_Metagenomics_2026/06_ANVIO/ERR5000343/METABAT_SUMMARY)
summary_folder_1=
summary_folder_2=

# copy fasta files
mkdir ALL_BINS
cp $summary_folder_1/bin_by_bin/*/*-contigs.fa ALL_BINS
cp $summary_folder_2/bin_by_bin/*/*-contigs.fa ALL_BINS

# run gtdbtk
gtdbtk \
  classify_wf \
  --genome_dir ALL_BINS \
  --out_dir GTDB \
  --extension fa \
  --cpus $SLURM_CPUS_PER_TASK \
  --pplacer_cpus 2
