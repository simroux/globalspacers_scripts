# How-to - reconstruct panels of Fig. 2
## Notes
Panel E is a schematic drawn separately. Data underlying panel A are further detailed in Fig. S7 and S8.  
## By panels:
### Panel A
Prepare the input file with `panel_a_prep.pl`  
Then Spacer_set_size_info.tsv can be used as input to prepare the panel A plot in R, as shown in `prepare_plots.R`  
### Panel B  
Prepare the input file with `prepare_alphadiv_plot_subsample.pl`  
Then spacer_sets_alphadiv_forplot_subsamples.tsv can be used as input to prepare the panel B plot in R, as shown in `prepare_plots.R`  
### Panel C
Prepare the input file compiling the ratios of singletons with `summarize_singleton_ratio.pl`  
Then input_figure_singletons.tsv can be used as input to prepare the panel C plot in R, as shown in `prepare_plots.R`  
### Panel D
Prepare the additional metadata file about each spacer set from the DuckDB database with `get_set_metadata_for_betadiv.sh`  
Use `gather_and_summarize_betadiv.pl` to process the different files, randomly subsample the sets and spacers within sets, and create the input file for the figure  
Finally, spacer_sets_betadiv_forplot_subsamples-summarized.tsv is used to prepare the top plots of panel D in R, as shown in `prepare_plots.R`  
In addition, `gather_and_summarize_betadiv_individuals.pl` can be used to process the files generated from groups of samples with the same individual or location information, and similarly randomly subsample the sets and spacers within sets to create the input file for the figure.  spacer_sets_betadiv_forplot_within-individuals.tsv is then used to prepare the bottom plots of panel D in R, as shown in `prepare_plots.R`  

