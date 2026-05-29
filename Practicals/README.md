# MBDP metagenomics – Practicals

__Table of Contents:__

1. [Introduction](#introduction)
2. [Setup](#setup)
3. [Data](#data)
4. [Quality control and trimming](#quality-control-and-trimming)
5. [Metagenome assembly](#metagenome-assembly)
6. [Assembly QC](#assembly-qc)
7. [Read-based taxonomy](#read-based-taxonomy)
8. [Viromics](#viromics)
9. [MAG binning in anvio](#mag-binning-in-anvio)
10. [Phylogenetic and functional analyses of MAGs](#phylogenetic-and-functional-analyses-of-mags)

## Introduction

During the course we will analyse metagenomic data from tundra soils collected in Kilpisjärvi, Finland.  
The whole dataset can be found in [ENA](https://www.ebi.ac.uk/ena/browser/view/PRJEB41762), and you can read the publication [here](https://doi.org/10.1186/s40793-022-00424-2).  
The samples were sequenced with both short-read (Illumina) and long-read (Nanopore) sequencing technologies, and for training purposes, we will focus on a small subset of of the data:  

- 2 Nanopore libraries (ERR5000342 and ERR5000343)
- 12 Illumina libraries (ERR4998593, ERR4998611, ERR4998615, ERR4998632, ERR4998657, ERR4998663, ERR4998600, ERR4998601, ERR4998602, ERR4998637, ERR4998638 and ERR4998640)

The Nanopore samples were also sequenced with Illumina, these are the matching accessions:  

| Sample | Vegetation | Nanopore   | Illumina   |
| -------|------------|------------|------------|
| m11216 | heathland  | ERR5000342 | ERR4998593 |
| m12208 | fen        | ERR5000343 | ERR4998600 |

## Setup

First create your own folder under the course project directory in Puhti:  

```bash
mkdir /scratch/project_2001499/$USER
```

And then clone this repository to your own folder:

```bash
cd /scratch/project_2001499/$USER
git clone https://github.com/MBDP-bioinformatics-courses/MBDP_Metagenomics_2026.git
```

## Data

Now create a folder for the course data and copy the data from the course project directory to your own directory:  

```bash
cd /scratch/project_2001499/$USER/MBDP_Metagenomics_2026
mkdir 01_DATA
cp -r /scratch/project_2001499/Data/* 01_DATA
```

This should copy all the data from the course project directory to your own directory. Note that the data is quite big, so it may take some time. After the data has been copied, check that both long and short reads were copied, and also the metadata file.  

## Quality control and trimming

### QC

We will start by checking the quality of the raw reads.  
For short reads, we will use FastQC and MultiQC, and for long reads NanoPlot.  

Before you start, allocate the computing resources with `sinteractive -i`.  
You will need 4 CPUs and 5 GB of memory, and it should not take more than 2 hours.  

__Short reads:__  
Run this once, it will analyse all the fastq files in the folder with FastQC and then summarize the results with MultiQC.  

```bash
mkdir 01_DATA/FASTQC

module load biokit
fastqc 01_DATA/Illumina/*.fastq.gz -o 01_DATA/FASTQC --threads $SLURM_CPUS_PER_TASK
module purge

module load multiqc
multiqc --interactive 01_DATA/FASTQC -o 01_DATA/FASTQC
module purge
```

__Long reads:__  
Run this command separately for both samples.  
Make sure to change the path to the fastq file and the name of the output folder.  
NanoPlot will give a warning about not finding Chrome, but it will still run. You can ignore the warning.  

```bash
/projappl/project_2001499/nano_tools/bin/NanoPlot \
    --threads $SLURM_CPUS_PER_TASK \
    -o path-to-output-folder \
    --only-report \
    --format png \
    --fastq path-to-nanopore-reads.fastq.gz
```

After QC is done, we will explore the results together.  

### Trimming

After checking the quality of the reads, we will remove the adapters from the short reads with `cutadapt`.  
The Nanopore reads will not be trimmed.  

```bash
module load cutadapt/4.9
mkdir 02_TRIMMED

for sample in 01_DATA/Illumina/*.R1.fastq.gz; do
    sample_name=$(basename $sample .novaseq.R1.fastq.gz)
    
    cutadapt \
        01_DATA/Illumina/${sample_name}.novaseq.R1.fastq.gz \
        01_DATA/Illumina/${sample_name}.novaseq.R2.fastq.gz \
        -o 02_TRIMMED/${sample_name}_R1_trimmed.fastq.gz \
        -p 02_TRIMMED/${sample_name}_R2_trimmed.fastq.gz \
        -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
        -A CTGTCTCTTATACACATCTGACGCTGCCGACGA \
        --minimum-length 50 \
        --cores $SLURM_CPUS_PER_TASK &> 02_TRIMMED/${sample_name}_cutadapt.log
done
```

After the trimiming is done, it would be good practice to check the quality of the trimmed reads again with FastQC and MultiQC.  
In case there is time, you can run both steps to the files in the `02_TRIMMED` folder and check the results.  

## Metagenome assembly

We will use three different approaches for metagenome assembly:  

- short-read assembly with `MEGAHIT`  
- long-read assembly with `Flye`  
- hybrid assembly with `metaspades`  

We will assemble only the samples were we have both short- and long-read data. So not all six short reads datasets.  
The assemblies will take some time, so you can prepare separate batch job scripts for each assembly approach and assemble always both samples in the same script. You can check the [CSC Puhti manual](https://docs.csc.fi/computing/running/creating-job-scripts-puhti/) or [CSC tutorial](https://csc-training.github.io/csc-env-eff/hands-on/batch_jobs/serial.html) on how to write a batch job script.  
The commands for each of the assemblies are given below. Check the options you used from the manual of each tool.  

```bash
mkdir 03_ASSEMBLY
```

```bash
/projappl/project_2001499/flye/bin/flye \
    --meta \
    --nano-raw 01_DATA/Nanopore/ERR5000342.nanopore.fastq.gz \
    --out-dir 03_ASSEMBLY/ERR5000342_flye \
    --threads $SLURM_CPUS_PER_TASK

/projappl/project_2001499/flye/bin/flye \
    --meta \
    --nano-raw 01_DATA/Nanopore/ERR5000343.nanopore.fastq.gz \
    --out-dir 03_ASSEMBLY/ERR5000343_flye \
    --threads $SLURM_CPUS_PER_TASK
```

```bash
module load spades/4.2.0

metaspades.py \
    -1 02_TRIMMED/ERR4998593_R1_trimmed.fastq.gz \
    -2 02_TRIMMED/ERR4998593_R2_trimmed.fastq.gz \
    --nanopore 01_DATA/Nanopore/ERR5000342.nanopore.fastq.gz \
    -o 03_ASSEMBLY/ERR5000342_hybrid \
    --threads $SLURM_CPUS_PER_TASK \
    --only-assembler

metaspades.py \
    -1 02_TRIMMED/ERR4998600_R1_trimmed.fastq.gz \
    -2 02_TRIMMED/ERR4998600_R2_trimmed.fastq.gz \
    --nanopore 01_DATA/Nanopore/ERR5000343.nanopore.fastq.gz \
    -o 03_ASSEMBLY/ERR5000343_hybrid \
    --threads $SLURM_CPUS_PER_TASK \
    --only-assembler
```

```bash
module load megahit

megahit \
    -1 02_TRIMMED/ERR4998593_R1_trimmed.fastq.gz \
    -2 02_TRIMMED/ERR4998593_R2_trimmed.fastq.gz \
    -o 03_ASSEMBLY/ERR5000342_megahit \
    --kmin 27 \
    --k-step 10 \
    --kmin-1pass \
    -t $SLURM_CPUS_PER_TASK

megahit \
    -1 02_TRIMMED/ERR4998600_R1_trimmed.fastq.gz \
    -2 02_TRIMMED/ERR4998600_R2_trimmed.fastq.gz \
    -o 03_ASSEMBLY/ERR5000343_megahit \
    --kmin 27 \
    --k-step 10 \
    --kmin-1pass \
    -t $SLURM_CPUS_PER_TASK
```

## Assembly QC

When the assemblies are ready, we will assess the assemblies with metaquast and choose the best approach for the downstream analyses.  

```bash
module load quast/5.2.0

metaquast.py \
    03_*/*_flye/assembly.fasta \
    03_*/*_hybrid/contigs.fasta \
    03_*/*_megahit/final.contigs.fa \
    -o 03_ASSEMBLY/QUAST \
    --max-ref-number 0
```

## Read-based taxonomy

Make a directory for read-based taxonomy & enter

```bash
mkdir /scratch/project_2001499/$USER/MBDP_Metagenomics_2026/04_TAXONOMY

cd /scratch/project_2001499/$USER/MBDP_Metagenomics_2026/04_TAXONOMY
```

Load the Metaphlan module & run Metaphlan using the array script after making any adjustments to the script if needed.

```bash
#!/bin/bash
#SBATCH --job-name=metaphlan
#SBATCH --account=project_2001499
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --array=0-11
#SBATCH --output=logs/metaphlan_%A_%a.out
#SBATCH --error=logs/metaphlan_%A_%a.err

module load metaphlan

THREADS=8
DB=/scratch/project_2001499/DBs/metaphlan
INPUT_DIR=/scratch/project_2001499/Data/Illumina/
OUT_DIR=./metaphlan_results

mkdir -p ${OUT_DIR}
mkdir -p logs

# Create array of R1 files
R1_FILES=(${INPUT_DIR}/*.R1.fastq.gz)

# Select current sample based on SLURM array task ID
R1=${R1_FILES[$SLURM_ARRAY_TASK_ID]}

# Extract sample name
SAMPLE=$(basename ${R1} .R1.fastq.gz)

# Define matching R2
R2=${INPUT_DIR}/${SAMPLE}.R2.fastq.gz

echo "Processing ${SAMPLE} ..."
echo "R1: ${R1}"
echo "R2: ${R2}"

metaphlan \
    ${R1},${R2} \
    --input_type fastq \
    --nproc ${THREADS} \
    --mapout ${OUT_DIR}/${SAMPLE}.mapout.txt \
    --db_dir ${DB} \
    -o ${OUT_DIR}/${SAMPLE}_profile.txt

echo "Finished ${SAMPLE}"
```

Merge files.

```bash
module load metaphlan
merge_metaphlan_tables.py ./metaphlan_results/*_profile.txt > merged_metaphlan.txt
```

Copy metadata to your own folder.

```bash
cp /scratch/project_2001499/Data/metadata.tsv .
```

Open an interactive session on Puhti with RStudio for 4h using the default settings; alternatively, use the small queue with 4 MB of memory, 4 cores, no NVMe, and 4h time.

Load the mia, miaViz and ggplot2 packages and set your working directory.

```r
library(mia)
library(miaViz)
library(ggplot2)
setwd("/scratch/project_2001499/$USER/MBDP_Metagenomics_2026/04_TAXONOMY")
```

1) Read from [OMA](https://microbiome.github.io/OMA/docs/devel/pages/import.html) and the import command's [help](https://microbiome.github.io/mia/reference/importMetaPhlAn.html) how to import Metaphlan objects.
Import data into an object called tse.

```r
tse <- mia::importMetaPhlAn("merged_metaphlan.txt", colData = sample_meta)
```

2) Inspect the treeSummarizedExperiment (TSE) object

```r
tse
assay(tse, "metaphlan")[1:3, 1:3]


rowData(tse) |> head()
# Check coldata, or sample data slot
colData(tse)
```

3) Check taxonomy ranks and how many unique phyla you have.

```r
getTaxonomyRanks()
getUnique(tse, rank = "phylum") |> head()
```

4) Let's visually check the abundance of the strains.

```r
plotAbundanceDensity(
    tse,
    layout = "jitter",
    assay.type = "metaphlan",
    n = 40, point.size = 1, point.shape = 19,
    point.alpha = 0.1
) +
    scale_x_log10(label = scales::percent)
```

5) Get top phylum and visualize.

```r
# Getting top taxa on a Phylum level
tse <- agglomerateByRank(tse, rank = "phylum")
top_taxa <- getTop(tse, top = 15, assay.type = "metaphlan")

# Inspect the top taxa
top_taxa

# Renaming the "phylum" rank to keep only top taxa and assign the rest to "Other"
phylum_renamed <- lapply(rowData(tse)$phylum, function(x) {
    if (x %in% top_taxa) {
        x
    } else {
        "Other"
    }
})
rowData(tse)$Phylum_sub <- as.character(phylum_renamed)

# Agglomerate the data based on specified taxa
tse_sub <- agglomerateByVariable(tse, by = "rows", f = "Phylum_sub")

# Visualizing the composition barplot, with samples ordered by the most abundant phylum
plotAbundance(
    tse_sub,
    assay.type = "metaphlan",
    order.row.by = "abund", order.col.by = "p__Pseudomonadota"
)
```

Alpha diversity

Read from https://microbiome.github.io/OMA/docs/devel/pages/alpha_diversity.html about the different alpha diversity indices.
Which one would you choose for this study?
Calculate all in one go using mia

```r
# The 'index' parameter allows computing multiple diversity indices
# simultaneously. Without specification, four standard indices are calculated:
# dbp_dominance, faith_diversity, observed_richness, and shannon_diversity.
tse <- mia::addAlpha(
    tse,
    assay.type = "metaphlan",
    detection = 10
)

```

06) Check alpha diversity by the vegetation type. If you have extra time, you can check how the numeric sample data correlates with Shannon, as exemplified here by moisture percentage. Which metric seems to have the highest correlation? You can check other alpha-diversity indices, too.

```r
library(patchwork)
library(scater)

# Create the plots
indices <- c(
    "dbp_dominance", "shannon_diversity"
)
plots <- lapply(
    indices,
    plotColData,
    object = tse,
    x = "vegetation",
    colour_by = "vegetation"
)

# Fine-tune visual appearance
plots <- lapply(
    plots, "+",
    theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()
    )
)

# Plot the figures
wrap_plots(plots, ncol = 1) +
    plot_layout(guides = "collect")


plotColData(tse, x = "shannon_diversity", y = "mositure_percent") +
    labs(x = "Shannon index", y = "Moisture %") +
    geom_smooth(method = "lm")

```

7) Check if the results are statistically significant (p<0.1) using linear models. pH here as example.


```r
library(dplyr)
df <- colData(tse) %>% as.data.frame()
lm(shannon_diversity ~ df$pH , data =df)   %>% summary()
```

8) Beta-diversity

Read about beta-diversity[https://microbiome.github.io/OMA/docs/devel/pages/community_similarity.html]

Which beta-diversity metric would you choose for this study?
Let's check Bray-Curtis and do unsupervised ordination analysis.

```r

# Run PCoA on the relabundance assay with Bray-Curtis distances
library(mia)
tse <- addMDS(
    tse,
    FUN = getDissimilarity,
    method = "bray",
    assay.type = "metaphlan",
    name = "MDS_bray"
)
```

9) Plot and color by vegetation. You can also color by the numeric variables. Which variable seems to drive dissimilarity between samples most?

```r
# Create ggplot object
p <- plotReducedDim(tse, "MDS_bray", colour_by = "vegetation")

# Calculate explained variance
e <- attr(reducedDim(tse, "MDS_bray"), "eig")
rel_eig <- e / sum(e[e > 0])

# Add explained variance for each axis
p <- p + labs(
    x = paste("PCoA 1 (", round(100 * rel_eig[[1]], 1), "%", ")", sep = ""),
    y = paste("PCoA 2 (", round(100 * rel_eig[[2]], 1), "%", ")", sep = "")
)

p
```

10) Let's do supervised ordination analysis with Bray-Curtis again using RDA.

```r
tse <- addRDA(
    tse,
    assay.type = "metaphlan",
    formula = assay ~ vegetation + pH + mositure_percent,
    distance = "bray",
    na.action = na.exclude
)
# Store results of PERMANOVA test

rda_info <- attr(reducedDim(tse, "RDA"), "significance")
```

11) Plot coloring by vegetation and add pH as a covariate.

```r
# Load packages for plotting function
library(miaViz)

# Generate RDA plot colored by clinical status
plotRDA(tse, "RDA", colour.by = "vegetation")
```

## Viromics

Make a directory for all virus analyses in your own directory (if not created yet):

```bash
cd /scratch/project_2001499/$USER/MBDP_Metagenomics_2026
mkdir 05_VIROMICS
```

NOTE: change ```$USER``` to your directory name.

All the following viromics work will be done within the ```05_VIROMICS``` directory.

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
mkdir GENOMAD
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
cd GENOMAD

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
cd ..
mkdir CHECKV
```

Run CheckV interactively:

```bash
cd CHECKV

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
cd ..
mkdir vOTUs
cd vOTUs

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

```bash
cd ..
mkdir IPHOP
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

### Annotating prokaryotic vOTUs with Phold

We will use [Phold](https://academic.oup.com/nar/article/54/1/gkaf1448/8415830), a protein structure-informed bacteriophage genome annotation tool for annotating the obtained vOTUs. (!) In a real project, you should use Phold only for prokaryotic viruses. Here, we keep all vOTUs at this step for learning purposes.

Note that Phold needs its database, which is already downloaded to ```/scratch/project_2001499/DBs/``` (and it's also specified in the batch job script below).

Make a directory ```PHOLD```for Phold output in your ```05_VIROMICS``` directory.

Make a batch job script (you can get a sample script from the updated repo, move to your ```src```directory and run ``` wget https://raw.githubusercontent.com/MBDP-bioinformatics-courses/MBDP_Metagenomics_2026/main/src/phold.sh``` or just copy this one below):

```bash
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
```
Check the used options from the manual.

Submit the batch job script. This run should take about 1 h.

**Phold output**

Explore the output:

- What files were created?
- Which PHROG categories are reported in *phold_all_cds_functions.tsv* ?
- What is PHROG?
- Did eukaryotic vOTUs also get phage-like predictions?

### IMG/VR v4 database

[IMG/VR v4](https://img.jgi.doe.gov/cgi-bin/vr/main.cgi) database: let's explore it online together! Separate instructions. 

Next version of the IMG/VR v4 db = [MetaVR](https://www.meta-virome.org/) database

Many more databases exist!

### Viromics summary

Let's think together what could be done in a real project with the obtained data. What other types of analyses could be run? Viral sequences in Kilpisjärvi soil samples were analysed in [Demina et al 2025](https://link.springer.com/article/10.1186/s40168-025-02053-6).

More about soil viruses:

- [A global atlas of soil viruses](https://www.nature.com/articles/s41564-024-01686-x) presents a comprehensive dataset compiled from almost 3K previously sequenced soil metagenomes -> about 38.5 K vOTUs.

- [Beneath the surface: Unsolved questions in soil virus ecology](https://www.sciencedirect.com/science/article/pii/S0038071725000732)

- [Soil viral diversity, ecology and climate change](https://www.nature.com/articles/s41579-022-00811-z)

If interested in other ecosystems, e.g. human gut, check [A genomic atlas of the human gut virome](https://www.biorxiv.org/content/10.1101/2025.11.01.686033v1).

**Other useful resources for future**

***Tools***

Many more tools for virus identification, annotation, and host prediction exist! See [Awesome-Virome](https://github.com/shandley/awesome-virome).

[Modular Viromics Pipeline](https://gitlab.com/ccoclet/mvp): nicely combines geNomad, CheckV, read mapping, functional annotation into a pipeline with several modules.

Be careful with AMGs 😊: [some guidelines](https://peerj.com/articles/11447/?utm_source=researchgate.net&utm_medium=article) and a [call for caution](https://www.nature.com/articles/s41564-025-02095-4). Upcoming: [CheckAMG](https://github.com/AnantharamanLab/CheckAMG) (pipeline under development).

***Webinars, meetings, conferences***

[European Virus Bioinformatics Center](https://evbc.uni-jena.de/) -> you can subscribe for newsletter, check annual ViBioM meetings, and a collection of virus bioinformatics tools

[ECR Viromics Webinar Series](https://coms.osu.edu/webinars/ecr-viromics-webinar-series), online, sign up to follow

[RNA Virus Journal Club](https://rdrp.io/journal-club/), online, sign up to follow, you can also nominate a speaker or even act as a chair!

[International Soil Virus Conference 2026](https://soilmicrobes.fr/international-soil-virus-conference-2026/): virtual participation may be still possible (?), 16-18 Jun 2026, France

[JGI VEGA symposium (Viral EcoGenomics and Applications)](https://jgi.doe.gov/work-with-us/events/vega-symposium), 18-19 Nov 2026, USA  

## MAG binning in anvio

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

```bash
cd /scratch/project_2001499/$USER/MBDP_Metagenomics_2026
mkdir 06_ANVIO
cd 06_ANVIO
```

In Puhti, we load our `anvi'o` installation with:  

```bash
export PATH="/projappl/project_2001499/anvio-dev/bin/:/projappl/project_2001499/anvio-gh/:$PATH"
export PYTHONPATH="/projappl/project_2001499/anvio-gh/:$PYTHONPATH"
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
You should check their online documentation, for example here:  https://anvio.org/help/main/programs/anvi-gen-contigs-database/.  
And since you're at it, familiarise yourself with two of the main `anvi'o` artifacts:  

- the `CONTIGS.db`: https://anvio.org/help/main/artifacts/contigs-db  
- the `PROFILE.db`: https://anvio.org/help/main/artifacts/profile-db  

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
  samtools sort sample.raw.bam
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
You can find an example script in `/scratch/project_2001499/$USER/MBDP_Metagenomics_2026/src/anvio.sh`.  
But remember:  

- You will need to modify the script to point the commands to your own files  
- You will need to run the script once for each assembly  

After you have submitted the script with `sbatch` the job will take a couple of hours to conclude.  
But once it is finished we are ready to bin the MAGs!  

### Creating and adding SSH keys

In this part we will need to connect to Puhti using a terminal emulator and the SSH protocol, the browser-based Puhti interface won't work.  
First, **if you haven't done this before**, you will need to create SSH keys and add them to my.csc.fi .  
Windows users will need to install the program `PuTTY` from the software centre, and Linux/Mac users will use the terminal app.  

- Windows: https://csc-training.github.io/csc-env-eff/hands-on/connecting/ssh-keys.html#windows
- Linux/Mac: https://csc-training.github.io/csc-env-eff/hands-on/connecting/ssh-keys.html#linux-and-macos

When your SSH keys have been added to my.csc.fi, you can now connect to Puhti:  

- Windows: https://docs.csc.fi/computing/connecting/ssh-windows
- Linux/Mac: https://docs.csc.fi/computing/connecting/ssh-unix

### Binning the MAGs with `anvi-interactive`  

For practical reasons, we will do this part using the program `screen`.  
This program acts as a virtual screen that keeps running when we logout of Puhti.  
You can read more about `screen` for example [here](https://www.geeksforgeeks.org/linux-unix/screen-command-in-linux-with-examples/), but what we need to know now is:  

```bash
# to launch a new screen
screen -S NAME_OF_THE_SCREEN

# to detach the screen hold the ctrl key 
# and then press a and d

# to see the list of the screens that are running
screen -ls

# and to retach a previous screen
screen -r NAME_OF_THE_SCREEN

# check if you are indeed inside a screen
# if the output is empty, you are NOT in a screen
echo $STY

# and to check in which login node we are
# to be able to retach the screen later
hostname -a
```

So, make sure that you are inside a screen and then:  

```bash
# connect interactively to a computing node
sinteractive -t 8:00:00 -m 24G -A project_2001499

# get the computing node id
# write down the output of this command, 
# we will need it for later
echo ${SLURMD_NODENAME}.bullx

# load the anvi'o installation
export PATH="/projappl/project_2001499/anvio-dev/bin/:/projappl/project_2001499/anvio-gh/:$PATH"
export PYTHONPATH="/projappl/project_2001499/anvio-gh/:$PYTHONPATH"

# get your port number
# also write down the output for later
grep $USER /scratch/project_2001499/anvio_ports.txt | cut -f 2

# now we launch anvi-interactive
# remember to adjust the command first
anvi-interactive \
  -c path-to-the-CONTIGS.db \
  -p path-to-the-merged-PROFILE.db \
  -P your-port-number \
  --ip-address 0.0.0.0 \
  --server-only
```

`anvi-interactive` is now hosting a web server that we can access using the browser (Chrome- or Firefox-based, preferentially).  
But to access this web server we need to do a small procedure called SSH tunelling.  
First, detach the screen with `ctrl + a + d` and logout of Puhti.  
Then, connect to Puhti again using a modified approach:  

**Linux/Mac**

```bash
ssh -L \
  PORT-NUMBER:COMPUTING-NODE-ID:PORT-NUMBER \
  YOUR-USER@LOGIN-NODE.csc.fi

# where:
# PORT-NUMBER: the output of 'grep $USER /scratch/project_2001499/anvio_ports.txt | cut -f 2'
# COMPUTING-NODE-ID: the output of 'echo ${SLURMD_NODENAME}.bullx'
# YOUR-USER: the output of 'echo $USER'
# LOGIN-NODE: the output of 'hostname -a'
```

Now you can retach the screen with `screen -r NAME_OF_THE_SCREEN`.  

**Windows with PuTTY**  

- Open PuTTY
- 'Host name' will now be the output of `echo ${SLURMD_NODENAME}.bullx`
- In the left menu go to Connection > SSH > Tunnels
  - 'Source port': the output of `grep $USER /scratch/project_2001499/anvio_ports.txt | cut -f 2`   
  - 'Destination': the output of `echo ${SLURMD_NODENAME}.bullx`, 
    followed by a colon (`:`),
    followed by the source port again,
    e.g. `r07c05.bullx:8123`
  - Click 'Add' and then 'Open'
- Now you can retach the screen with `screen -r NAME_OF_THE_SCREEN`.  

**Binning the MAGs (finally!)** 



## Phylogenetic and functional analyses of MAGs

