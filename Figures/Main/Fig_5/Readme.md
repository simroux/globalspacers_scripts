# How-to - reconstruct panels of Fig. 5
## Notes
Data underlying panels A and B are further detailed in Fig. S15 and Fig. S16.
## By panels:
### Panels A and B
Panels A and B are based on the analysis of virus sequences from previous studies in which possible broad targeting (i.e. targeting by phylogenetically distinct hosts) was observed. Some of these virus sequences are not in IMG/VR, so we need to identify and post-process spacer hits to these viruses, as illustrated in "Analyses/Target_Custom_datasets/".  
Next, we prepare the input files for each panel as follows:
```
./fig_prep_cover_for_network.pl -p a
./fig_prep_cover_for_network.pl -p b
```
This will create input files for edges and nodes for each panel in folders named "panel_a" and "panel_b", that can then be imported in Cytoscape for visualization. 
### Panel C
The right-side schematic were drawn manually. To obtain the numbers indicated to the left of the schematic:  
We first extract the corresponding data from the database. Warning: this requires a lot of memory !  
```
## First, we link each relevant spacer to the corresponding lca genus and lca origin of this array if of the right confidence, ignoring cases with no genus
CREATE TABLE tmp_table AS (
    SELECT cluster_id, spacer_filt_clusters.spacer_id, crispr_array, lca_origin, lca_genus FROM spacer_filt_clusters JOIN 
        (SELECT spacer_id, crispr_array, lca_origin, lca_genus FROM spacer_filt_tbl JOIN 
            (SELECT repeat_cluster, lca_origin, lca_genus FROM array_tbl WHERE lca_genus!='NA' AND (lca_origin='Genome_medium-confidence' OR lca_origin='Genome_high-confidence')) AS selected_arrays
            ON (selected_arrays.repeat_cluster=spacer_filt_tbl.crispr_array)) AS selected_spacers
        ON (spacer_filt_clusters.spacer_id=selected_spacers.spacer_id)
    )

####
## Then we join this with the spacer hits information, summarise informoation for each virus - taxon pairs, and export
COPY(
    SELECT * FROM 
        (SELECT target_id, lca_genus, first(lca_origin) as lca_origin, count(distinct(crispr_array)) as n_repeats, count(distinct(tmp_table.cluster_id)) as n_clusters FROM tmp_table JOIN imgvr_hits_filt ON (tmp_table.cluster_id=imgvr_hits_filt.cluster_id) GROUP BY target_id, lca_genus)
    WHERE n_clusters >= 10)
TO 'Target_to_taxon.tsv' (HEADER, DELIMITER '\t');

```  

Next, run `./select_relevant_taxa_pairs.pl run` to summarize the data for each virus, and the code in `count_multitaxa.R` to get the different numbers. 
Finally, `./get_genus_pairs_for_multiclass.pl` can be used to identify pairs of genera that tend to (co-)target the same viruses.  

