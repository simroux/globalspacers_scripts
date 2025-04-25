# How-to - matching spacers to custom virus datasets and extract spacer hits and repeat information
**Note: This uses functions available in [SpacerExtractor](https://code.jgi.doe.gov/SRoux/spacerextractor), version 0.9.5 was used in this study.**  

## Matching spacers to new virus sequences
**Note: here we use the viruses presented in Fig. 5, including some that are part of IMG/VR, but we reprocess as this also provides a good illustration of how such an analysis can be conducted on any virus dataset.**.

XX LINK TO NOTEBOOK ? XX

First, we build the required database for the potential targets, here "selected_uvigs.fna".  
```
conda activate spacerextractor
spacerextractor create_target_db -i ../../Data/Additional_data/selected_uvigs_for_network.fna -d selected_uvigs_db -t 4 --replace_spaces
```

Next, we map all spacers to these viruses, using the file "nr_spacers_hq.fna" we generated for the IMG/VR and IMG/PR analysis:
```
spacerextractor map_to_target -i ../Target_IMGVR_IMGPR/nr_spacers_hq.fna -d selected_uvigs_db -o All_spacers_vs_selected_viruses/ -t 4
```

## Post-process spacer matches and select relevant hits
XX Python to get information XX

