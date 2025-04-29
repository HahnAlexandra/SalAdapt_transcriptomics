#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw629/transcriptomics/log-files
#SBATCH --job-name=salmon_quant
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=48:00:00
#SBATCH --output=%x.%A.out
#SBATCH --error=%x.%A.err
#SBATCH --partition=base

cd $WORK/transcriptomics/data/X208SC24034727-Z01-F001/03.Trimming

for READ in $(ls | cut -f 1-3 -d "_" | sort | uniq); do

  salmon quant -i $WORK/transcriptomics/tonsa_transcriptome/salmon_index \
	  -l A \
	  -1 ${READ}_1.trim.fq.gz -2 ${READ}_2.trim.fq.gz \
	  --seqBias \
	  --gcBias \
	  -p 16 \
	  -o ../04.Alignment/${READ}_quant

done
