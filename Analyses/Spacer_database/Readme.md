# Documentation for scripts in 'Analyses/Spacer_Database' folder
## Database building
`create_duckdb.sh`
This includes the commands, run into duckdb, that import data about samples, repeats, and spacer, from the tsv files into a duckdb database.  
These input data include:
* Runs_to_ecosystem_and_sequencing_and_study_for_db-Jul28-24.tsv: Information about SRA entries used in the study (see "Analyses/Extra/SRA_metadata")
* Array_info_filtered_for_db-Oct24-25.tsv: Information about individual repeats (note: 'array' may be used for 'repeat' in some places throughout, sorry..) (see "Analyses/Repeat_database")
* All_spacers_info_filtered-Jul19-24.tsv: Information about spacers extracted from SRA data. Obtained from SpacerExtractor (see [https://code.jgi.doe.gov/SRoux/spacerextractor](https://code.jgi.doe.gov/SRoux/spacerextractor))
* All_spacers_info_filtered_clusters-Jul19-24.tsv: Information about the global clustering (100% identity over 100% of the length) of spacers across the entire study, performed with "custom_clustering_spacers.py"

`custom_clustering_spacers.py`
Python script used to cluster spacers based on exact match over their entire length

## Diversity of spacers within set (Fig. 2, Figs. S7, S8, S9)
To evaluate the diversity of spacers observed within spacer sets, the steps include:
* Export relevant data from duckdb, as seen in `extract_alphadiv_data.sh`
* Calculate diversity for each set with `calculate_alpha_diversity.pl`
Note: to avoid requiring too much memory, this computation is split into two batches, one with sets of 10 spacers or more, and the other with sets including less than 10 spacers. Both files are extracted separately from duckdb (see `extract_alphadiv_data.sh`) and need to be processed separately afterwards (`calculate_alpha_diversity.pl -m`). Finally, the two result files must be concatenated as follows:
```
cp spacer_sets_alphadiv_up10.tsv spacer_sets_alphadiv.tsv
cat spacer_sets_alphadiv_md10.tsv | grep -v "n_spacers" >> spacer_sets_alphadiv.tsv
```

Clustering of spacer sets a 95% and 80% identity can be done as in `cluster_by_sample.pl`.

## Identification of common spacers
This is useful later, and more convenient to compute once and use the list later (since most spacesr are not common)  
`./identifiy_common_spacers.pl run`  

## Linking spacers to ecosystem
This is also useful later and more convenient to compute once here. This is in two steps, first an overall count of ecosystem type is generated for each unique spacer from duckdb (see `extract_alphadiv_data.sh`). Then, the script `link_spacer_to_eco_final.pl` can be used to clean up the (rare but existing) complicated cases where a single spacer is found in multiple ecosystems.  
