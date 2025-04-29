#!/bin/bash
#SBATCH --job-name=FastP
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=60G
#SBATCH --time=02:00:00
#SBATCH --output=fastp.out
#SBATCH --error=fastp.err

cd $WORK/transcriptomics/data/X208SC24034727-Z01-F001/01.RenamedData

for NAME in $(ls | cut -f 1-6 -d "_" | sort | uniq); do

  BASENAME=$(echo $NAME | cut -f 2,3,6 -d "_")

  fastp -w 10 -Q -g -i ${NAME}_1.fq.gz -I ${NAME}_2.fq.gz -o ../03.Trimming/${BASENAME}_1.trim.fq.gz -O ../03.Trimming/${BASENAME}_2.trim.fq.gz

done
