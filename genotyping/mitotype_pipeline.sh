#!/bin/bash

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
   echo "pipeline to determine mitotypes of A tonsa using muscle for alignment and mr bayes."
   echo
   echo "Syntax: ./mitotype_pipeline.sh input_fasta plot_name"
   echo
   echo "Mandatory input:"
   echo "input_fast         file containing the new samples to analyze. In fasta format. This file should be in the main directory with all scripts"
   echo "plot_name          name for output plot for this run"
   echo
   echo "This script requires a number of fasta files containing the reference sequences from figuerra et al in addition to the scripts in the main directory. Do not modify these files"
   echo
   echo "all output will be directed to ./output"
  exit 0
fi

if [[ ! -f "$1" ]]; then
    echo "The fasta file "$1" does not exist. Exiting."
    exit
fi

if [[ ! -e output ]]; then
    mkdir output
fi

#we will start with fasta files. these need to be aligned, then run mr bayes

# raw fasta is: /data/copepods/mitotypes/reference_haplotypes.fasta


#merge new samples with reference samples
cat "$1" \
    reference_haplotypes.fasta \
    > ./output/to_align.fasta


# align all data
#conda activate muscle
echo "alignment starting"
muscle -align ./output/to_align.fasta \
                                        -output ./output/aligned.fasta
echo "alignment done"
echo " "

# make mr bayes file
## need to have ape and stringr installed
echo "making mr bayes file"

Rscript make_mr_bayes.R

# run mr bayes

# https://www.ebi.ac.uk/Tools/services/rest/muscle/result/muscle-I20200720-155944-0426-80668140-p1m/aln-fasta
echo " "
echo "Mr bayes starting. This will take 10 or 15 minutes"

/opt/homebrew/opt/mrbayes/bin/mb ./output/tonsa_mb.nex > ./output/mr_bayes_log.txt

if grep -q "Error" ./output/mr_bayes_log.txt; then
    echo "ERROR in Mr Bayes! Are there spaces in your fasta names? \nExiting."
    exit
fi


echo " "
echo "Mr bayes done\n"
echo " "
# HKY + I + G for sub model
# for HKY: nst=2
# I + G = invariable site plus discrete Gamma model: rates=invgamma


#plot output
echo "generating plot"
Rscript plot_tree.R

mv ./output/tree_plot.pdf ./output/"$2"_plot.pdf



