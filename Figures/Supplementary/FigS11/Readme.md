# How-to - reconstruct panels of Fig. S11
First, we extract the relevant viruses and spacer hits for taxa of interests  
The list of taxa of interest is provided in "list_taxa.tsv"  
DuckDB commands in `extract_spacer_from_taxa.sh` will generate 3 files:  
- selected_viruses.tsv: with information about IMG/VR v4 sequences selected based on our taxa of interest
- selected_spacer_hits.tsv: with information about spacer hits to these viruses
- selected_spacer_info.tsv: with information about spacers involved in these spacer hits

Next, `select_hit_info.pl` can be used to generate the file "fig_hits_input_vr-atypical_phages.tsv" that will provide an overview of hits on these viruses of interest, and `link_hits_to_taxa-ecosystem.pl` can be used to generate the file "Input_for_atypical_taxa_ecosystem.tsv" that will provide information on the taxa and ecosystem of the hits to atypical phages and archaeoviruses.  

Finally, the different panels can be prepared in R using the code included in `prepare_plots.R`.  

