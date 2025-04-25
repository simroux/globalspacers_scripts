# How-to - reconstruct panels of Fig. 4
## Notes
Some of the input files must be generated with scripts included in "Analyses/Target_IMGVR_IMGPR/Target_coverage/" and "Analyses/Target_IMGVR_IMGPR/Beyond_near_exact/".
## By panels:
### Panel A
Input file must be generated via `get_distribution_spacernumber.pl`. Next, the panel can be prepared in R using `plot_coverage.R`.  
### Panel B
Input file is generated with `get_sample_coverage_summary.pl`. Next, the panel can be prepared in R using `plot_coverage.R`.  
### Panels C & D
Input files are generated in the folder "Analyses/Target_IMGVR_IMGPR/Beyond_near_exact/". Next, the plot in panels C and D can be prepared in R as indicated in `plot_coverage.R`.  
### Panel E & G
Input file is generated with `count_targeting_type_vs_virus_features.pl`. Next, the plots in panel E and G can be prepared in R as indicated in `plot_coverage.R`.  
### Panel F
Input file is generated with `prepare_known_host_figure.pl`. Next, the plot in panel F can be prepared in R as indicated in `plot_coverage.R`.  