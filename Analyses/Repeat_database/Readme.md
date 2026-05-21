# Documentation for scripts in 'Analyses/Repeat_Database' folder
## Preparation of repeat sequence database (for documentation purposes)
This folder includes scripts used to prepare the repeat sequence database. The file "", provided in "Spacerdb_raw_files.tar.gz" (see main Readme), is already the end result of these different steps, so the scripts are provided only for documentation purposes. The differents steps include:
- `./link_repeat_to_taxon_and_cluster.pl`: first clustering of the predicted repeats, selection of a representative, and association with genomes (and genome taxonomy) when available.
- `python check_repeat_complexity.py`: read a fasta file of clustered repeats, and flag sequences of low complexity
- `./filter_restriction_and_quality.pl`: process the full list of repeats, adjust LCA confidence when needed, add LCA information from metagenome contigs, and output an updated list
- `./reassign_crispr_type.pl`: re-run repeatType on the selected repeats, and update the predicted CRISPR type in the list.
- `./link_cas_prediction_to_repeat.pl`: process the CRISPRCasTyper results to link individual repeats to neighboring Cas loci, will generate the file "cas_by_repeat_cluster"
- `./reassign_crispr_type_with_cas.pl`: clean up and integrate the Cas-based CRISPR type predictions to the repeat cluster table, will generate "Array_info_filtered_for_db-updated_type_and_cas.tsv"


The final file generated corresponds to "Array_info_filtered_for_db-Apr23-26.tsv

