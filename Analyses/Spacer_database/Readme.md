# Documenation on scripts available in 'Analyses/Spacer_Database' folder
## Database building
`create_duckdb.sh`
This includes the commands, run into duckdb, that import data about samples, repeats, and spacer, from the tsv files into a duckdb database.  
These input data include:
* Runs_to_ecosystem_and_sequencing_and_study_for_db-Jul28-24.tsv: Information about SRA entries used in the study (see "Analyses/Extra/SRA_metadata")
* Array_info_filtered_for_db-Nov1-24.tsv: Information about individual repeats (note: 'array' may be used for 'repeat' in some places throughout, sorry..) (see "Analyses/Repeat_database")
* All_spacers_info_filtered-Jul19-24.tsv: Information about spacers extracted from SRA data. Obtained from SpacerExtractor (see [https://code.jgi.doe.gov/SRoux/spacerextractor](https://code.jgi.doe.gov/SRoux/spacerextractor))
* All_spacers_info_filtered_clusters-Jul19-24.tsv: Information about the global clustering (100% identity over 100% of the length) of spacers across the entire study, performed with "custom_clustering_spacers.py"

`custom_clustering_spacers.py`
Python script used to cluster spacers based on exact match over their entire length

