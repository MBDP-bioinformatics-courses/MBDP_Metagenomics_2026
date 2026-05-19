#!/bin/bash
#SBATCH --job-name=anvio
#SBATCH --account=project_2001499
#SBATCH --time=4:00:00
#SBATCH --mem=24G
#SBATCH --partition=small
#SBATCH --cpus-per-task=20
#SBATCH --gres=nvme:48
#SBATCH --output=%x-%j.out

# initialise
module load anvio

# you need to add these below:
# assembly_id: name of your assembly (e.g. ERR5000342)
# output_dir: path to the folder where the results will be stored (e.g. /scratch/project_2001499/$USER/06_ANVIO/ERR5000342)
# contigs_fasta: path to the fasta file with the assembled contigs (e.g. /scratch/project_2001499/$USER/02_ASSEMBLY/ERR5000342_flye/assembly.fasta)
# illumina_data: path to the folder containing the trimmed Illumina data (e.g. /scratch/project_2001499/$USER/02_TRIMMED)
assembly_id=''
output_dir=''
contigs_fasta=''
illumina_data=''

# create output folder
mkdir $output_dir && cd $_

# reformat fasta
anvi-script-reformat-fasta $contigs_fasta \
                           -o CONTIGS.fa \
                           -r CONTIGS-reformat.txt \
                           -l 5000 \
                           --prefix $assembly_id \
                           --simplify-names

# create contigs db
anvi-gen-contigs-database -f CONTIGS.fa \
                          -o CONTIGS.db \
                          -n $assembly_id \
                          -T $SLURM_CPUS_PER_TASK

# rum hmms for SSUs and SCGs
anvi-run-hmms -c CONTIGS.db \
              -T $SLURM_CPUS_PER_TASK

# get SCG taxonomy
anvi-run-scg-taxonomy -c CONTIGS.db \
                      -T $SLURM_CPUS_PER_TASK

# map the illumina reads

## create bowtie index
bowtie2-build CONTIGS.fa CONTIGS.idx --threads $SLURM_CPUS_PER_TASK &> bowtie2-build.log

## map the samples
for r1 in ${illumina_data}/*_R1_trimmed.fastq.gz
do
  sample=`basename $r1 _R1_trimmed.fastq.gz`
  r2=${illumina_data}/${sample}_R2_trimmed.fastq.gz

  bowtie2 -1 $r1 \
          -2 $r2 \
          -S ${sample}.sam \
          -x CONTIGS.idx \
          -p $SLURM_CPUS_PER_TASK \
          --no-unal

  samtools view -F 4 -bS -@ $SLURM_CPUS_PER_TASK ${sample}.sam -o ${sample}-RAW.bam
  samtools sort -@ $SLURM_CPUS_PER_TASK ${sample}-RAW.bam -o ${sample}.bam 
  samtools index -@ $SLURM_CPUS_PER_TASK ${sample}.bam
  rm ${sample}.sam ${sample}-RAW.bam
done

# create profile dbs
for bamfile in *.bam
do
  sample=`basename $bamfile .bam`

  anvi-profile -i $bamfile \
               -c CONTIGS.db \
               -o ${sample}-PROFILE \
               -S $sample \
               -T $SLURM_CPUS_PER_TASK \
               --skip-hierarchical-clustering
done

# merge profiles
anvi-merge *-PROFILE/PROFILE.db \
           -c CONTIGS.db \
           -o MERGED \
           -S $assembly_id \
           --enforce-hierarchical-clustering
