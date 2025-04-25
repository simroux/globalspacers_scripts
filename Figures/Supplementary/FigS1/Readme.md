# Benchmarking consistency of taxonomic assignment of repeats across genomes/bins
Note: this uses some information available in "Data/Additional_data/", and you will need to unzip the tsv file the first time around.
## Input preparation for top panel
`./get_lca_rank.pl run`  
Compile the statistics on the LCA rank for repeats, based only on genome taxonomy (isolate and bins)  
## Input preparation for bottom panel
`./test_error_ngenomes.pl run`  
Run the subsampling and estimate the consistency of taxonomy for individual repeats for different number of genomes considered.   
## Figures
`plot_fig_lca_repeat.R`
in R, to generate the figure  