#!/bin/bash
#SBATCH --job-name=FastQC
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G
#SBATCH --time=00:30:00
#SBATCH --output=$HOME/transcriptomics/log-files/fastqc.out
#SBATCH --error=$HOME/transcriptomics/log-files/fastqc.err

cd $WORK/transcriptomics/data/X208SC24034727-Z01-F001/01.RawData

find . -name "*.fq.gz" -print0 | xargs -0 -P 10 fastqc -t 10 -o ../02.QC
