#!/bin/bash
#SBATCH --job-name=XX
#SBATCH --output=%x-%j.out
#SBATCH --account=project_2001499
#SBATCH --time=XX:00:00
#SBATCH --mem=XXG
#SBATCH --partition=small
#SBATCH --cpus-per-task=16

## Remember to change each XX to something meaningful above

# Load necessary modules
module load YOUR_MODULE/XX

# Run your command here

