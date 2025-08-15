# pipeline to analyze mitotype data

This pipeline holds an automated method to analyze mitotype data for A. tonsa roughly based on the method from of Figueroa et al., 2020: Phylogeography of Acartia tonsa Dana, 1849 (Calanoida: Copepoda)and phylogenetic reconstruction of the genus Acartia Dana, 1846.

`combined.fasta` contains the COI sequence data for Baltic and North Sea individuals used in the analysis.


---

View the help file with: `./mitotype_pipeline.sh --help` or `./mitotype_pipeline.sh -h`

To reproduce the pipeline as in the manuscript, run the following:

```bash
./mitotype_pipeline.sh combined.fasta tree_plot

```

This script requires a number of fasta files containing the reference sequences in addition to the scripts in the main directory. Do not modify these files. These files are:

- `reference_haplotypes.fasta`: the reference sequences of species to compare with
- `make_mr_bayes.R`: R script to make Mr Bayes input
- `plot_tree.R`: R script to generate final tree
 
