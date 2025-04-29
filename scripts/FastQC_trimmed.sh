#!/bin/bash
#SBATCH --job-name=FastQC
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G
#SBATCH --time=00:30:00
#SBATCH --output=fastqc_trimmed.out
#SBATCH --error=fastqc_trimmed.err

cd $WORK/transcriptomics/data/X208SC24034727-Z01-F001/03.Trimming

fastqc -t 10 -o . *.trim.fq.gz
