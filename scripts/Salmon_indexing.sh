#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw629/transcriptomics/log-files
#SBATCH --job-name=salmon_index
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --output=%x.%A.out
#SBATCH --error=%x.%A.err
#SBATCH --partition=base

salmon index -t $WORK/transcriptomics/tonsa_transcriptome/trinity.trimmomatic.above500.noPhiX.fasta.gz -i $WORK/transcriptomics/tonsa_transcriptome/salmon_index -k 31
