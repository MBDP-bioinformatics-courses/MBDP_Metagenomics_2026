# MBDP metagenomics – Practicals

__Table of Contents:__

1. [Introduction](#introduction)
2. [Setup](#setup)
3. [Data](#data)
4. [Quality control](#quality-control)
5. [Metagenome assembly](#metagenome-assembly)
6. [Assembly QC](#assembly-qc)
7. [Read-based taxonomy](#read-based-taxonomy)
8. [Viromics](#viromics)
9. [Genome-resolved metagenomics](#genome-resolved-metagenomics)
10. [MAG QC and taxonomy](#mag-qc-and-taxonomy)
11. [MAG annotation](#mag-annotation)
12. [Automatic binning](#automatic-binning)

## Introduction

During the course we will analyse metagenomic data from soils collected in Kilpisjärvi, Finland. The data originates from the publication: [ZZZ](LINK) and are available in SRA under accession number XXX.  
The samples were sequenced with both short-read (Illumina) and long-read (Nanopore) sequencing technologies, but for training purposes, we will focus on only a subset of of the data, 2 nanopore sequenced samples (ERR5000342 and ERR5000343) and 6 Illumina sequenced samples (ERR5000342, ERR5000343, ERR5000344, ERR5000345, ERR5000346, and ERR5000347). The data are already available on Puhti for the course analyses.
The matching samples sequenced with both technologies are ...

## Setup

First create your own folder under the course project directory on Puhti, e.g.:

```bash
mkdir /scratch/project_2001499/$USER
```

NOTE: change ```$USER``` to your directory name.

After that, clone this github reposoitory to your own folder.  

```bash
cd /scratch/project_2001499/$USER   
git clone https://github.com/MBDP-bioinformatics-courses/MBDP_Metagenomics_2026.git
```

## Data

When the course github pages are cloned, create folder for the course data and copy the data from the course project directory to your own directory:
We will make separate folders for short- and long-read data.

```bash
mkdir /scratch/project_2001499/$USER/MBDP_Metagenomics_2026/01_DATA/short_read
mkdir /scratch/project_2001499/$USER/MBDP_Metagenomics_2026/01_DATA/long_read

cp /scratch/project_2001499/Data/Illumina/* /scratch/project_2001499/$USER/MBDP_Metagenomics_2026/01_DATA/short_read/
cp /scratch/project_2001499/Data/Nanopore/* /scratch/project_2001499/$USER/MBDP_Metagenomics_2026/01_DATA/long_read/
```

Copy also the metadata file to your own directory:

```bash
cp /scratch/project_2001499/Data/metadata.tsv /scratch/project_2001499/$USER/MBDP_Metagenomics_2026/01_DATA/
```

After copying, verify that you have the right files in your data folders with ```ls```.  

## Quality control and trimming

### QC

We will start by checking the quality of the raw reads. For short reads, we will use FastQC and MultiQC, and for long reads, we will use NanoPlot and nanoQC.  

Before you start, allocate resources with `sinteractive -i`. You will need 4-6 CPUs and 2-5 GB of memory for the QC. This should not take more than 2 hours.

__Short reads:__  
Run this once, it will analyse all the fastq files in the folder with FastQC and then summarize the results with MultiQC.  

```bash
module load biokit
fastqc 01_DATA/short_read/*.fastq.gz -o 01_DATA/FASTQC --threads $SLURM_CPUS_PER_TASK
multiqc --interactive 01_DATA/FASTQC -o 01_DATA/
module purge
```

__Long reads:__  
Run these to commands with both samples separately, changing the path to the fastq file and output folder name.

```bash
/projappl/project_2001499/nano_tools/bin/NanoPlot \
  -o nanoplot_out -f png --fastq path-to/your_raw_nanopore_reads.fastq

/projappl/project_2001499/nano_tools/bin/nanoQC -o nanoQC_out path-to/your_raw_nanopore_reads.fastq
```

After QC is done, we will explore the results together.  

### Trimming

After checking the quality of the reads, we will remove adapter from short reads with cutadapt.  
The long reads will not be trimmed.  

```bash
module load cutadapt/4.9

for sample in 01_DATA/short_read/*.R1.fastq.gz; do
    sample_name=$(basename $sample .novaseq.R1.fastq.gz)
    cutadapt \
        -a FW_ADAPTER -A RV_ADAPTER \
        -o 02_TRIMMED/${sample_name}_R1_trimmed.fastq.gz \
        -p 02_TRIMMED/${sample_name}_R2_trimmed.fastq.gz \
        01_DATA/short_read/${sample_name}_R1.fastq.gz \
        01_DATA/short_read/${sample_name}_R2.fastq.gz \
        --cores $SLURM_CPUS_PER_TASK \
        --minimum-length 50 \
        > 00_LOGS/${sample_name}_cutadapt.log
done
```

## Metagenome assembly

We will use three different approaches for metagenome assembly: short-read assembly, long-read assembly, and hybrid assembly. For short-read assembly, we will use MEGAHIT, for long-read assembly, we will use Flye, and for hybrid assembly, we will use metaspades.  

We will assemble only the samples were we have both short- and long-read data. So not all six short reads datasets.  
The assemblies will take some time, so you can prepare separate batch job scripts for each assembly approach and assemble always both samples in the same script. You can check the [CSC Puhti manual](https://docs.csc.fi/computing/running/creating-job-scripts-puhti/) on how to write a batch job script.  
The commands for each of the assemblies are given below. Check the options you used from the manual of each tool.  

```bash
/projappl/project_2001499/flye/bin/flye \
    --meta \
    --nano-raw 01_DATA/long_read/ERR5000342.nanopore.fastq.gz \
    --out-dir 03_ASSEMBLY/ERR5000342_flye \
    --threads $SLURM_CPUS_PER_TASK

/projappl/project_2001499/flye/bin/flye \
    --meta \
    --nano-raw 01_DATA/long_read/ERR5000343.nanopore.fastq.gz \
    --out-dir 03_ASSEMBLY/ERR5000343_flye \
    --threads $SLURM_CPUS_PER_TASK
```

```bash
module load spades/4.2.0

metaspades.py \
    -1 02_TRIMMED/ERR4998593_R1_trimmed.fastq.gz \
    -2 02_TRIMMED/ERR4998593_R2_trimmed.fastq.gz \
    --nanopore 01_DATA/long_read/ERR5000342.nanopore.fastq.gz \
    --only-assembler -o 03_ASSEMBLY/ERR5000342_hybrid \
    --threads $SLURM_CPUS_PER_TASK

metaspades.py \
    -1 02_TRIMMED/ERR4998600_R1_trimmed.fastq.gz \
    -2 02_TRIMMED/ERR4998600_R2_trimmed.fastq.gz \
    --nanopore 01_DATA/long_read/ERR5000343.nanopore.fastq.gz \
    --only-assembler \
    -o 03_ASSEMBLY/ERR5000343_hybrid \
    --threads $SLURM_CPUS_PER_TASK
```

```bash
module load megahit

megahit \
    -1 02_TRIMMED/ERR4998593_R1_trimmed.fastq.gz \
    -2 02_TRIMMED/ERR4998593_R2_trimmed.fastq.gz \
    -o 03_ASSEMBLY/ERR5000342_megahit \
    --kmin 27 --k-step 10 --kmin-1pass \
    -t $SLURM_CPUS_PER_TASK

megahit \
    -1 02_TRIMMED/ERR4998600_R1_trimmed.fastq.gz \
    -2 02_TRIMMED/ERR4998600_R2_trimmed.fastq.gz \
    -o 03_ASSEMBLY/ERR5000343_megahit \
    --kmin 27 --k-step 10 --kmin-1pass \
    -t $SLURM_CPUS_PER_TASK
```

## Assembly QC

When the assemblies are ready, we will assess the assemblies with metaquast and choose the best approach for the downstream analyses.  

```bash
module load quast/5.2.0

metaquast.py \
    03_ASSEMBLY/*_flye/assembly.fasta \
    03_ASSEMBLY/*_hybrid/contigs.fasta \
    03_ASSEMBLY/*_megahit/final.contigs.fa \
    -o 03_ASSEMBLY/QUAST \
    --max-ref-number 0
```

## Read-based taxonomy

## Viromics

Make a directory for all virus analyses in your own directory (if not created yet):

```bash
mkdir /scratch/project_2001499/$USER/MBDP_Metagenomics_2026/05_VIROMICS
```

NOTE: change ```$USER``` to your directory name.

### Identifying viral contigs using geNomad

There are many different tools for predicting viral contigs from metagenomes. In this course, we will use geNomad. [Read about it](https://www.nature.com/articles/s41587-023-01953-y) and check its [GitHub pages](https://github.com/apcamargo/genomad). Good documentation also [here](https://portal.nersc.gov/genomad/pipeline.html). How does it work?

Note that geNomad needs its database, which is already downloaded to ```/scratch/project_2001499/DBs/``` (and it's also specified in the batch job script below).

**Running geNomad**

In your 05_VIROMICS directory, create a sample list (*sample_list.txt*), which will have sample names:

```bash
ERR5000342
ERR5000343
```

Make a directory for geNomad output:

```bash
mkdir /scratch/project_2001499/$USER/MBDP_Metagenomics_2026/05_VIROMICS/GENOMAD/
```

You can find a sample batch job script (*genomad.sh*) in ```/scratch/project_2001499/$USER/MBDP_Metagenomics_2026/src/```, but check all paths and change if needed:

```bash
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
```

Submit the job:

```bash
sbatch /path-to/genomad.sh /path-to/sample_list.txt
```

Check the used options (and others) by ```genomad -h``` and ```genomad end-to-end -h```. Remember to enter ```export PATH="/projappl/project_2001499/genomad/bin:$PATH"```first.

Note that we used a while loop in this example, how does it work? What other options do you have if you need to run the same tool for multiple samples?

**geNomad output**

Explore the output you got. Have a look at the log files first:  

- What steps did geNomad run?  
- How many viral contigs were identified in each sample? How about plasmids?  

Find summary tables for each sample, where viral contigs are listed:  

- What viral taxa were predicted?  
- Any RNA viruses? Can they be here?
- What length do viral contigs have? What length is a bacteriophage genome on average vs other viral groups such as giant viruses?  
- Are there proviruses predicted?  

### Quality control with CheckV

We will assess the quality and completeness of viral contigs identified by geNomad with CheckV. [Read about the tool](https://www.nature.com/articles/s41587-020-00774-7) and [how it works](https://bitbucket.org/berkeleylab/checkv/src/master/).

Note that CheckV needs its database, which is already downloaded to ```/scratch/project_2001499/DBs/``` (and it's also specified in the batch job script below).

Before running CheckV, we can combine geNomad viral contigs (fna files) from two samples into one set. Since some contigs may have same names in both samples, we should add a sample-based prefix first to contig names so that all headings are unique in a combined fna file:

```bash
cd /scratch/project_2001499/$USER/MBDP_Metagenomics_2026/05_VIROMICS/GENOMAD/

sed "s/^>/>ERR5000342_/" ERR5000342/assembly_summary/assembly_virus.fna > ERR5000342_virus.fna

sed "s/^>/>ERR5000343_/" ERR5000343/assembly_summary/assembly_virus.fna > ERR5000343_virus.fna

cat *_virus.fna > virus_combined.fna
```

Check that prefixes were added with e.g. ```head```command and you can also check that your combined file contains the right number of sequences with ```seqkit stats```:

```bash
module load biokit

seqkit stats virus_combined.fna
```

**Running CheckV**

Make a directory for CheckV analyses:

```bash
mkdir /scratch/project_2001499/$USER/MBDP_Metagenomics_2026/05_VIROMICS/CHECKV/
```

Run CheckV interactively:

```bash
cd /scratch/project_2001499/$USER/MBDP_Metagenomics_2026/05_VIROMICS/CHECKV/

sinteractive -A project_2001499 -m 10G -c 8 

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

export PATH="/projappl/project_2001499/checkv/bin:$PATH" 

checkv end_to_end \
/scratch/project_2001499/$USER/MBDP_Metagenomics_2026/05_VIROMICS/GENOMAD/virus_combined.fna \
virus_combined_checkv.out \
-d /scratch/project_2001499/DBs/checkv-db-v1.5/ \
-t $SLURM_CPUS_PER_TASK
```

Check the options you used with ```checkv -h``` and ```checkv end_to_end -h```.

Running CheckV will take a few minutes.

**CheckV output**

Explore the output, especially the summary file *quality._summary.tsv*:

- Are there any proviruses predicted? Are those contigs that were flagged as proviruses by geNomad now also flagged as proviruses by CheckV? Any new proviruses compared to geNomad predictions?*
- Are most contigs of low, medium or high quality? See how the quality corresponds to completeness.
- Are there any 100% complete viral genomes listed?
- Are there contigs with kmer_freq > 1? This indicates that the viral genome is represented multiple times in the contig, which is quite rare.
- Any warnings?

*Note that we used the geNomad output file as the input for CheckV. geNomad output .fna has proviral sequences already cut from host-derived flanking regiongs. These contigs are also renamed with provirus coordinates (e.g. ERR5000343_contig_14802|provirus_18530_57545). CheckV may not recognize them as proviruses, since they don't have host genes anymore. In some cases, however, it does flag them as proviral and thus suggests cutting them more. In addition, it can also mark as proviruses some contigs that were not recognized as proviral by geNomad. CheckV also outputs the proviruses.fna file, which contains proviral sequences cut as CheckV suggests (note contig names changes, e.g. ERR5000342_contig_11875_1 766-4342/4342). This means that in a real project, one would need to run CheckV again on proviruses for getting correct data on proviral_length, gene_count, viral_genes, host_genes, etc. Note that in the second round, CheckV may "want" to cut even more from some proviruses, but usually, it's not run more than twice.  

 (!) In a real project, CheckV output is typically used for filtering some predictions out. Common thresholds for metagenomic viral contigs include:  

- at least 1 virus gene identified by CheckV;
- host to virus gene count ratio no more than 1:1;
- length minimum of 5 kbp or 10 kbp, unless a genome is >=50% complete (but not shorter than 1 kbp anyway).
  

Different thresholds are used for metatranscriptomes.

In this course, we won't filter any viral contigs.

### Dereplicating viral contigs into vOTUs

Since some viral contigs could have been present (and assembled) in both samples, the datasets from the two samples may overalp. To dereplicate viral contigs into viral operational taxonomic units = vOTUs, which roughly correspond to viral species, we can use BLAST (as [parallel BLAST at CSC](https://docs.csc.fi/apps/blast/#usage-of-pb-parallel-blast-at-csc)) and anicalc.py and aniclust.py scripts from CheckV. Standard thresholds for dereplicating into vOTUs: 95% average nucleotide identity and 85% alignment fraction. Check [Minimum Information about an Uncultivated Virus Genome (MIUViG)](https://www.nature.com/articles/nbt.4306) for more info on how vOTU is defined.

For training purposes, we can use all viral contigs predicted by geNomad (without addittional CheckV-based filtering) for dereplicating as follows:

```bash
# make a directory for vOTUs

mkdir /scratch/project_2001499/$USER/MBDP_Metagenomics_2026/05_VIROMICS/vOTUs

cd /scratch/project_2001499/$USER/MBDP_Metagenomics_2026/05_VIROMICS/vOTUs

module load biokit 

# blast viral contigs against themselves

pb blastn -dbnuc ../GENOMAD/virus_combined.fna -query ../GENOMAD/virus_combined.fna \
-outfmt '6 std qlen slen' -max_target_seqs 1000 -out virus_combined.tsv

# pb blast will take about 25 min

module load biopythontools

# calculate ANI values

python /projappl/project_2001499/anicalc.py -i virus_combined.tsv -o virus_combined_ani.tsv

# cluster contigs 

python /projappl/project_2001499/aniclust.py --fna ../GENOMAD/virus_combined.fna \
--ani virus_combined_ani.tsv --out virus_combined_clusters.tsv \
--min_ani 95 --min_tcov 85 --min_qcov 0

# save the first column of the tsv with clusters into a txt file => vOTUs IDs

cut -f1 virus_combined_clusters.tsv > vOTUs_IDs.txt

# extract vOTU sequences based on their IDs from the original fasta file

seqtk subseq ../GENOMAD/virus_combined.fna vOTUs_IDs.txt > vOTUs.fna

# check how the final vOTUs fasta files looks like

seqkit stats vOTUs.fna
```

How many viral contigs were predicted from two samples in total? How many were retained as vOTUs?

### Linking vOTUs to putative hosts

We will use iPHoP for linking vOTUs to their putative bacterial and archaeal hosts (note: not suitable for eukaryotic viruses). Although some viral contigs were classified as eukaryotic viruses in our dataset, we'll still include them here for training puprposes, but in a real project, you should exclude them from the iPHoP input.

iPHoP integrates multiple methods for host predictions: which ones? [Check the publication](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3002083) and [documentation](https://bitbucket.org/srouxjgi/iphop/src/main/). What are these methods based on?

Note that iPHoP needs its database, which is already downloaded to ```/scratch/project_2001499/DBs/``` (and it's also specified in the batch job script below). We will run the latest version of the default database, but it is also possible to construct your own database by adding e.g. MAGs obtained from the same samples (see "Adding bacterial and/or archaeal MAGs to the host database" in [documentation](https://bitbucket.org/srouxjgi/iphop/src/main/)).

**Running iPHoP**

Make a directory for iPHoP output:

```
mkdir /scratch/project_2001499/$USER/MBDP_Metagenomics_2026/05_VIROMICS/IPHOP
```

Sample batch job script (found in ```/scratch/project_2001499/$USER/MBDP_Metagenomics_2026/src/```):

```bash
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
```
Submit with ```sbatch```.

Check the used options from manual (how to call it?).

**iPHoP output**

Explore the output. Find the *Host_prediction_to_genus_m75.csv* file. Note that we have applied the 75% cut-off threshold, which is OK for family-level predictions, but the 90% cut-off threshold should be used for genus-level predictions. In a real project, you would need to filter the predictions based on these thresholds.

- How many vOTUs got predictions (% of total)? Are there multiple predictions for some vOTUs?
- Which methods are listed as used ones?
- How many valid genus- vs family-level predictions are there? How about higher levels?
- How do host predictions match read-based taxonomic profiles of the samples? I.e., are most abundant bacterial/archaeal taxa among the predicted hosts?

### Further reading

Let's think together what could be done in a real project with the obtained data. What other types of analyses could be run? Viral sequences in Kilpisjärvi soil samples were analysed in [Demina et al 2025](https://link.springer.com/article/10.1186/s40168-025-02053-6).

More about soil viruses:

- [A global atlas of soil viruses](https://www.nature.com/articles/s41564-024-01686-x) presents a comprehensive dataset compiled from almost 3K previously sequenced soil metagenomes -> about 38.5 K vOTUs.

- [Beneath the surface: Unsolved questions in soil virus ecology](https://www.sciencedirect.com/science/article/pii/S0038071725000732)

- [Soil viral diversity, ecology and climate change](https://www.nature.com/articles/s41579-022-00811-z)

If interested in other ecosystems, e.g. human gut, check [A genomic atlas of the human gut virome](https://www.biorxiv.org/content/10.1101/2025.11.01.686033v1), or for marine viruses, check the [Tara Oceans project](https://www.tara-oceans-science.org/viruses/).

### IMG/VR v4 database

[IMG/VR v4](https://img.jgi.doe.gov/cgi-bin/vr/main.cgi) database: let's explore it online together! Separate instructions. 

Next version of the IMG/VR v4 db = [MetaVR](https://www.meta-virome.org/) database

Many more databases exist! Also specific ones, like [PaVE](https://pave.niaid.nih.gov/) for papillomaviruses.

### Other useful resources for future

**Tools**

Many more tools for virus identification, annotation, and host prediction exist! See [Awesome-Virome](https://github.com/shandley/awesome-virome).

[Modular Viromics Pipeline](https://gitlab.com/ccoclet/mvp): nicely combines geNomad, CheckV, read mapping, functional annotation into a pipeline with several modules.

Be careful with AMGs 😊: [some guidelines](https://peerj.com/articles/11447/?utm_source=researchgate.net&utm_medium=article) and a [call for caution](https://www.nature.com/articles/s41564-025-02095-4). Upcoming: [CheckAMG](https://github.com/AnantharamanLab/CheckAMG) (pipeline under development).

**Webinars, meetings, conferences**

[European Virus Bioinformatics Center](https://evbc.uni-jena.de/) -> you can subscribe for newsletter, check annual ViBioM meetings, and a collection of virus bioinformatics tools

[ECR Viromics Webinar Series](https://coms.osu.edu/webinars/ecr-viromics-webinar-series), online, sign up to follow

[RNA Virus Journal Club](https://rdrp.io/journal-club/), online, sign up to follow, you can also nominate a speaker or even act as a chair!

[International Soil Virus Conference 2026](https://soilmicrobes.fr/international-soil-virus-conference-2026/): virtual participation may be still possible (?), 16-18 Jun 2026, France

[JGI VEGA symposium (Viral EcoGenomics and Applications)](https://jgi.doe.gov/work-with-us/events/vega-symposium), 18-19 Nov 2026, USA  

## Genome-resolved metagenomics

Microbial genomes allow us to study in detail things like microbial metabolism, structural variation, horizontal gene transfer, etc.  
Ideally, we would like to work with complete, circular genomes—**can you think of a reason for this?**  
**But wait, is this even achievable for metagenomic data?**  

Let's take a look again at the report from `MetaQuast`, which gave us information on the metagenome assemblies:  

- **How long is the longest contig in each assembly?**
- **What about the N50 values?**

Based on these, do you think that:  

1. **Each contig represents a complete bacterial or achaeal genome**; or 
2. **The genomes are most likely fragmented into many contigs, each covering only a fraction of the complete genome?**

If you have answered **1**, congratulations, you can skip this part!  

But if you have answered **2**, which is more likely the correct answer, you have then realised that long-read technologies do not guarantee complete, chromosome-level microbial genomes in a single contiguous sequence.  
It is in fact likely that our metagenome assemblies contain mostly genomes that are fragmented across many contigs.  
We can't really improve genome contiguity in metagenomes without much, much longer reads.  
**But could we somehow group the contigs that are coming from the same population to obtain a better representation of their genomes?**  

It turns out we can, and there are different bioinformatic strategies to achieve this.  
We will use mainly the concepts of **sequence similarity** and **differential coverage**:  

- two contigs that come from the same original population will share similar nucleotide composition  
- and their coverage signals will be similar across different samples

We can use these two simple concepts to identify the contigs that originate from the same population and group them into a genomic bin.  
We usually refer to these bins as "population genomes" or "metagenome-assembled genomes" (MAGs).  
We will mostly use `anvi'o`, which is an open-source, community-driven **an**alysis and **vi**sualization platform for microbial **'o**mics (https://anvio.org).  
Although `anvi'o` does include automatic binning programs (e.g. `MetaBat2`), we will focus on manual, interactive binning.  

Let's start by making a directory for the genome-resolved analyses:  

```
mkdir /scratch/project_2001499/$USER/06_ANVIO
```

In Puhti, you can load the `anvi'o` environment with:  

```bash
module load anvio
```

For each assembly, we need to run several commands to prepare the files for `anvi'o`.  
Let's go through them first **without running anything for now**; we will do this later using `sbatch`.  
The four commands below will take the fasta file containing the assembled contigs and process the sequences in many ways:  

```bash
anvi-script-reformat-fasta contigs.fasta
anvi-gen-contigs-database -f reformatted-contigs.fasta
anvi-run-hmms -c CONTIGS.db
anvi-run-scg-taxonomy -c CONTIGS.db
```

**What is each command doing?**  
You should check their online documentation, for example here:  https://anvio.org/help/8/programs/anvi-gen-contigs-database.  
And since you're at it, familiarise yourself with two of the main `anvi'o` artifacts:  

- the `CONTIGS.db`: https://anvio.org/help/8/artifacts/contigs-db  
- the `PROFILE.db`: https://anvio.org/help/8/artifacts/profile-db  

The four commands above will create the `CONTIGS.db` and populate it with information about the sequences, such as the location of open reading frames, SSU rRNA genes and single-copy genes that are used to assess genome quality.  

The `PROFILE.db` stores information about the contigs across multiple samples, including nucleotide coverage and variability.  
The commands below will loop through the Illumina reads and map them to the assembled contigs:  

```bash
# first we index the contigs
bowtie2-build reformatted-contigs.fasta

# then we map the samples
for sample in sample1 sample2 sample3
do
  bowtie2 -1 R1_reads.fastq -2 R2_reads.fastq
  samtools view sample.sam
  samtools sort sample.bam
  samtools index sample.bam
done
```

Finally,  to create the `PROFILE.db`:  

```bash
# first we loop through the bam files and make them into a PROFILE.db
for sample in sample1 sample2 sample3
do
  anvi-profile sample.bam
done

# then we merge them into a single PROFILE.db
anvi-merge PROFILE1.db PROFILE2.db PROFILE3.db
```

**What are these commands doing specifically?**  
**What is the difference between `anvi-profile` and `anvi-merge`?**

Once everything is well understood, you are ready to run these commands using `sbatch`.  
You can find an example script in `/scratch/project_2001499/$USER/MBDP_Metagenomics_2026/src/anvio_sbatch.sh`.  
But remember:  

- You will need to modify the script to point the commands to your own files  
- You will need to run the script once for each assembly  

After you have submitted the script with `sbatch` the job will take a couple of hours to conclude.  
But once it is finished we are ready to bin the MAGs!

## MAG QC and taxonomy

## MAG annotation

## Automatic binning
