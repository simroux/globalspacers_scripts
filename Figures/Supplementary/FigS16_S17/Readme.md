# How-to - reconstruct panels of Fig. S16 and Fig. S17
The scripts for figures S16 and S17 are combined as these two figures are based on the same visualization, and correspond to virus-taxon connections presented in panels C and D of Fig. 5, respectively.  
The data are based on analyses described in "Analyses/Target_Custom_datasets/", along with the input files generated for the network visualizations for Fig. 5 ("Figures/Main/Fig_5/"). The following commands will reconstruct the input files needed for the R plots:
```
./combine_all_hits_for_figs.pl -p c
./combine_all_hits_for_figs.pl -p d
```
with the first one generating the input file for Fig. S16, and the second the input file for Fig. S17.  
Next, the commands listed in `plot_spacer_coverage.R` can be used to generate the coverage plots, which will be exported separately for each virus.
