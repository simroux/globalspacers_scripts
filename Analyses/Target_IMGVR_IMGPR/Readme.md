# How-to - matching spacers against IMG/VR and IMG/PR
**Note: this can now be done (much) more easily via the functions available in [SpacerExtractor](https://code.jgi.doe.gov/SRoux/spacerextractor), and is only provided here as documentation of the exact scripts that were used in the study.**  

## Matching spacers to IMG/VR and IMG/PR 
First, we need to build bowtie-1 compatible databases from IMG/VR and IMG/PR. IMG/VR and IMG/PR nucleotide fasta files can be downloaded from the [IMG website](https://img.jgi.doe.gov).  
`./build_db.pl -i IMG_VR_sequences-high_confidence.fna -d Database_IMGVR/` (same command with adjusted inputs and outputs for IMG/PR).

Next, we can export all spacers from the DuckDB database and create a non-redundant fasta file of spacers  
`export_all_spacers.sh`

Finally, we can use bowtie1 to map spacers to these potential targets, with a Perl script wrapper to process the whole database and post-process the sam file in order to get a nice tsv result file in the end (**again, this is now better handled in [SpacerExtractor](https://code.jgi.doe.gov/SRoux/spacerextractor), and is only provided here for documentation purposes**).  
`./spacermap.pl -i nr_spacers_hq.fna -d Database_IMGVR/ -o ./results_nr_imgvrhq/ -t /tmp/tmp_searches/ -c 8 -q`
The same approach can be used to match spacers to IMGPR by adjusting the input and output

## Post-process spacer matches and select relevant hits
Next, we post-process the output of spacermap to remove spurious hits and prepare the tsv file to be imported in the duckdb database.  
First, we look for low-complexity regions in the database(s) with dustmasker.  
`dustmasker -in IMG_VR_sequences-high_confidence.fna -out Database_IMGVR/IMG_VR_sequences-high_confidence.fna.dustmasker`

Next, because we work with such large database (IMG/PR and IMG/PR) we predict CRISPR arrays in these genomes and mask the hits from these loci. Matching spacer to its own array is not necessarily an issue per se, but in our case since we wanted to report the number of "spacer-protospacer" hits, we opted to remove them. Below are the commands for IMG/VR, and the same can be done for IMG/PR. These are based on a combination of de-novo detection (minced) and reference-based detection (matching to repeats in our repeat database).
```
cd Arrays_in_VR_PR/
minced -gff ../IMG_VR_sequences-high_confidence.fna IMG_VR_sequences-high_confidence.minced.txt
./repeatmap.pl -i /.../SE_Db_r214GBIMG_0.9/Repeats.fna -d ../Database_IMGVR/ -o repeats-vs-imgvr/
./summarize_regions_to_mask.pl -i IMG_VR_sequences-high_confidence.minced.txt -m repeats-vs-imgvr/Repeats_vs_IMGVR_db_20_all.regions.tsv -o ../Database_IMGVR/IMG_VR_sequences-high_confidence.fna.predicted_arrays.tsv
```

Now we are ready to post-process the spacer hits we had found in the first step, and try to remove as many spurious hits as we can.  
`./post_process.pl -d Database_IMGVR/ -o results_nr_imgvrhq/nr_spacers_hq_vs_IMGVR_db_all.matches.tsv.gz`

Finally, we can gather all these hits (including if they are split across multiple files, e.g. because we processed the spacers by separate batches), and change the column headers to be easier to work with in duckdb (also removing some columns that are not useful)  
`./gather_all_best_matches.pl -i results_nr_imgvrhq/nr_spacers_hq_vs_IMGVR_db_all.matches_clean.tsv.gz -c ../../Data/Spacer_db/All_spacers_info_filtered_clusters-Jul19-24.tsv -o Results_nr_imgvrhq_all_best_matches.tsv`

`All_spacers_info_filtered_clusters-Jul19-24.tsv` is generated based on the overall spacer dataset, as described in `Analyses/Spacer_database/`.
This entire process is applied to both IMG/VR and IMG/PR, and leads to the generation of `Results_nr_imgvrhq_all_best_matches-Sept13-24.tsv` and `Results_nr_imgprhq_all_best_matches-Sept13-24.tsv`, both available in the data bundle (`Data/').

## Import hits into the duckdb database
The steps used to import the hits to IMG/VR and IMG/PR, along with the information about IMG/VR and IMG/PR sequences, are listed in  
`import_to_duckdb.sh`

## Run PAM detection based on hits to IMG/VR and IMG/PR
The corresponding scripts are available in the folder "PAM_detection".  

First, we need to export the information for sequences upstream and downstream from spacer hits from DuckDB. Here for IMG/VR first, then for IMG/PR.  
```
COPY (SELECT crispr_array, upstream, downstream, COUNT(*) FROM (SELECT crispr_array, cluster_id FROM spacer_hits_imgvr GROUP BY crispr_array, cluster_id) AS tmp, imgvr_hits_filt WHERE imgvr_hits_filt.cluster_id=tmp.cluster_id GROUP BY crispr_array, upstream, downstream ORDER BY crispr_array) TO 'Neighborhood_for_pam_motif_all_imgvr.tsv' (HEADER, DELIMITER '\t');
COPY (SELECT crispr_array, upstream, downstream, COUNT(*) FROM (SELECT crispr_array, cluster_id FROM spacer_hits_imgpr GROUP BY crispr_array, cluster_id) AS tmp, imgpr_hits_filt WHERE imgpr_hits_filt.cluster_id=tmp.cluster_id GROUP BY crispr_array, upstream, downstream ORDER BY crispr_array) TO 'PAM_viruses_plasmids/Neighborhood_for_pam_motif_all_imgpr.tsv' (HEADER, DELIMITER '\t');
```

Next, we proceed with a de novo identification of conserved residues upstream and downstream of spacer hits 
`./detect_conserved_positions.pl -d ../../../Data/`  
Then, we evaluate whether these conserved residues are consistent with known PAMs (or variants of known PAMs), and assign each repeat to the most likely PAM (previously described or, when needed, newly identified)
`./match_conservation_to_known_pams.pl -d ../../../Data/`  
This will generate the files "Stat_motif_detection.tsv" and "Detailed_stats_motifs.tsv" used in Figure 3 and Supplementary Figure 12  


## Get stats on the number of hits by spacer
We can export the number of hits (i.e. number of distinct vOTUs or PTUs) hit by cluster, both for all sequences and for only high-quality UViGs or near-complete plasmids. This is done from the duckdb database, as described in `export_spacer_hit_stats.sh`  
We can get further stats on the ANI distance between high-quality genomes targeted by the same spacer, using the script `get_distrib_dist_by_cluster.pl`. Specifically, to generate ANI data used in Fig. 3E, we can use `./get_distrib_dist_by_cluster.pl -i ../../Data/Additional_data/selected_subset_spacers_fordist.tsv -c 8`  

## Get information on coverage of targets by spacers
The corresponding scripts are available in the folder "Target_coverage".  

First we need to extract the corresponding information from duckdb, using commands listed in `extract_spacer_info.sh`.  
Next, we can use the scripts `get_uvig_coverage.pl` and `get_plasmid_coverage.pl` to calculate the coverage of individual targets (virus or plasmids) by spacers linked to individual repeats. And use `add_sp_info_to_coverage.pl` to add some more information about the type of hits. Note that this last script must be run twice, once for viruses as `./add_sp_info_to_coverage.pl run` and once for plasmids as `./add_sp_info_to_coverage.pl run -p`.  

We can also use `get_uvig_coverage-bysample.pl` to get per-sample coverage for all pairs in which the overall coverage is at least 10% of the UViG or 1kb.

## Import and process additional hits (i.e. hits with 2 and 3 mismatches) for relevant target-repeat pairs
The corresponding scripts are available in the fodler "Beyond_near_exact". They rely on files generated in the "Target_coverage" folder, as well as the full results from `spacermap.pl`.  

First, we generate a tsv file for import into duckdb that includes all spacer hits at 2 and 3 mismatches for high-quality targets with at least one repeat shoring 10 spacer hits or 200bp target coverage, with the script `collect_extra_hits.pl`. Note that this file is also pre-computed and provided as part of the "Data" package.  

Next, we use the commands lised in `extract_relevant_extra_hits.sh` to create tsv files listing hits and spacer information for relevant virus-repeat pairs. This will generate the following files:  
```
All_hits_vr_hq_2-3mis.tsv
All_spacers_vr_hq_2-3mis.tsv
All_hits_pr_hq_2-3mis.tsv
All_spacers_pr_hq_2-3mis.tsv
```  

Next, we use the scripts `get_uvig_to_repeat_profiles.pl` to generate tsv files with information about the number of hits at different mismatch levels for virus-repeat pairs.  

## Evaluate hits for phages with known hosts
The corresponding scripts are available in the folder "Known_hosts".  

First, we use `export_hits.py` to extract information for relevant phages, using the list of phages with known hosts provided in "Data/Additional_data/". The command should be `python export_hits.py -u ../../../Data/Additional_data/Phages_with_known_hosts.tsv -d ../../Spacer_database/global_spacer.duckdb`.  

Next, we use `./cross_reference_known_hosts_and_CRISPR_hits.pl run` to compare the taxonomy of the targeting microbe to the one for the host listed in the database.

Finally, back in the "Beyond_near_exact" folder, we run `./add_virus_info_to_profiles.pl run` to add this information about known host to the targeting profile table.  

## Add coverage and mismatch information for viruses targeted by multiple taxa  
This is based on a list of viruses identified as targeted by repeats assigned to distinct taxa (classes or above). This list can be generated based on scripts provided as part of the main figure 5, but is also provided as "Data/Additional_data/List_multiclass_uvigs_nohq.txt".  

First, in "Target_coverage/", we export the spacer hits corresponding to these viruses from duckdb using the code indicated in `extract_spacer_info.sh`. Next, we calculate the coverage for each virus-repeat pair and add additional information as follows:
```
./get_uvig_coverage.pl -m run
./add_sp_info_to_coverage.pl -m run
```
The "-m" parameter will make the scripts work on these "multitaxa" uvigs, and lead to the generation of "uvigadditional_multitaxa_coverage_by_spacer.tsv" and "uvigadditional_multitaxa_coverage_by_spacer-with_sp_info.tsv". 

Next, in "Beyond_near_exact/", we can collect additional hits (i.e. 2 or 3 mismatches) for these viruses with `./collect_extra_hits.pl -m run`. This will generate the file "Additional_hits_imgvr-multitaxa-Dec17.tsv", already provided in the data bundle. We can then use this file to extract the corresponding spacers and hits from duckdb using the code provided in `extract_relevant_hits.sh`. This will generate the files "All_hits_vr_additional_multitaxa_2-3mis.tsv" and "All_spacers_vr_additional_multitaxa_2-3mis.tsv", that can then be used as input for `./get_multitaxa_uvig_to_repeat_profiles.pl run`. This will generate the file "Multitaxa_virus_to_array_hits_profile.with_cc.tsv", which includes information on mismatch profiles for these additional uvigs.  


