# How-to - reconstruct panels of Fig. 4
## Notes
Some of the input files must be generated with scripts included in "Analyses/Target_IMGVR_IMGPR/Target_coverage/" and "Analyses/Target_IMGVR_IMGPR/Beyond_near_exact/".
## By panels:
### Panel A
Input file must be generated via `get_distribution_spacernumber.pl`. Next, the panel can be prepared in R using `plot_coverage.R`.  
### Panel B
Input file is generated with `prepare_known_host_figure.pl`. Next, the plot can be prepared in R as indicated in `plot_coverage.R`.  
### Panel C
Input file is generated in the folder "../../../Analyses/Target_IMGVR_IMGPR/Known_hosts/. Next, the plot can be prepared in R as indicated in `plot_coverage.R`.  
### Panels D & E
Input file is generated in the folder "../../../Analyses/Target_IMGVR_IMGPR/Beyond_near_exact/. Next, the plots can be prepared in R as indicated in `plot_coverage.R`.  
