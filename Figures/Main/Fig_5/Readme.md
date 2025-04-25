# How-to - reconstruct panels of Fig. 1
## Notes
Some parts of panels A and B are schematics drawn separately. Data underlying panel C are further detailed in Fig. S15 and Fig. S16.
## Data pre-processing
We first export from duckdb the number of spacers linking individual virus targets with individual repeats, along with the repeat LCA, in the file "Target_to_repeat.tsv", using commands listed in `prep_export_multitaxa.sh`. Next, we identify viruses targeted by multiple taxa with the script `identify_multitaxa.pl`, which will generate the files "Detailed_multitaxa.tsv", "Summary_multitaxa.tsv", as well as "List_multiclass_uvigs.txt" and "List_multiclass_uvigs_nohq.txt".  
For both panels A and panels B, we need information about spacer hits coverage and mismatch profiles for these selected viruses. This information is already available for high-quality viruses, as it was generated as part of Fig. 4. The section "Add coverage and mismatch information for viruses targeted by multiple taxa" for the readme in "Analysis/Target_IMGVR_IMGPR" includes information about how to generate the required files for the other viruses detected as multitaxa but not already processed.  
Once these files have been generated, we can prepare a new summary that incorporate these coverage and mismatch profiles information with `./add_cover_and_mismatch_to_multitaxa_summary.pl run`, which will generate the file "Summary_multiclass_wprofile.tsv".
## By panels:
### Panel A
The top part is a schematic drawn separately. 
The input file "Multitaxa_frequency.tsv" for the bottom part can be generated with `./get_stats_panelA.pl run`, and the figure panel can then be created in R, based on the code in `plot_multitaxa.R`. This will also create a file called "List_of_interest.tsv" which is reused to prepare the data for panel B.  
### Panel B
The left part are schematics drawn separately.  
For the right part, the input file "Multitaxa_features.tsv" can be generated with `./get_stats_panelB.pl run`, and the figure is plotted in R using the code in `plot_multitaxa.R`.  
### Panels C and D
Panels C and D are based on the analysis of virus sequences from previous studies in which possible broad targeting (i.e. targeting by phylogenetically distinct hosts) was observed. Some of these virus sequences are not in IMG/VR, so we need to identify and post-process spacer hits to these viruses, as illustrated in "Analyses/Target_Custom_datasets/".  
Next
XXX WE NEED TO MAKE SCRIPTS TO PREPARE THE NETWORKS FOR EACH PANEL, BASED ON THE LIST OF VIRUSES IN EACH PANEL XXX

