# How-to - matching spacers to custom virus datasets and extract spacer hits and repeat information
**Note: This uses functions available in [SpacerExtractor](https://code.jgi.doe.gov/SRoux/spacerextractor), version 0.9.5 was used in this study.**  

## Matching spacers to new virus sequences
**Note: here we use the viruses presented in Fig. 5, including some that are part of IMG/VR, but we reprocess as this also provides a good illustration of how such an analysis can be conducted on any virus dataset.**.


First, we build the required database for the potential targets, here "selected_uvigs_for_network.fna".  
```
conda activate spacerextractor
spacerextractor create_target_db -i ../../Data/Additional_data/selected_uvigs_for_network.fna -d selected_uvigs_db -t 4 --replace_spaces
```

Next, we map all spacers to these viruses, using the file "nr_spacers_hq.fna" we generated for the IMG/VR and IMG/PR analysis:
```
spacerextractor map_to_target -i ../Target_IMGVR_IMGPR/nr_spacers_hq.fna -d selected_uvigs_db -o All_spacers_vs_selected_viruses/ -t 4
```

## Post-process spacer matches and select relevant hits
From there, we want to extract the taxonomy and other information about spacers with hits to our viruses of interest. This id done with the script `extract_spacers_from_hits.py`, as follows:  
```
python extract_spacers_from_hits.py --in_file All_spacers_vs_selected_viruses/nr_spacers_hq_vs_selected_uvigs_db_all_hits.tsv --out_file All_spacers_vs_selected_viruses/nr_spacers_hq_vs_selected_uvigs_db_all_hits_spacer_info.tsv -d ../Spacer_database/global_crispr_db.duckdb
```

**Note: A more complete examples of how to interact with the spacer database based on hits to a given target is available in the [example notebooks](XX XX TO FILL THE URL XX XX) presented with the SpacerDB release.**
