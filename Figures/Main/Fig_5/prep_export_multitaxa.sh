## For taxon <-> UViG analysis
# First we need to connect uvigs to individual taxa
# For that, we link each virus-spacer hit to the arrays that include this spacer, including the lca genus and lca origin of this array, ignoring cases with no genus
CREATE TABLE tmp_table AS (SELECT cluster_id, spacer_filt_clusters.spacer_id, crispr_array, lca_origin, lca_genus FROM spacer_filt_clusters JOIN (SELECT spacer_id, crispr_array, lca_origin, lca_genus FROM spacer_filt_tbl JOIN (SELECT repeat_cluster, lca_origin, lca_genus FROM array_tbl WHERE lca_genus!='NA') AS selected_arrays ON (selected_arrays.repeat_cluster=spacer_filt_tbl.crispr_array)) AS selected_spacers ON (spacer_filt_clusters.spacer_id=selected_spacers.spacer_id));
# Then we export the virus - array pairs while counting the number of distinct spacers, and we include info about the confidence of the taxo
COPY(SELECT target_id, crispr_array, first(lca_origin), first(lca_genus), COUNT(tmp_table.cluster_id) as n_clusters FROM tmp_table JOIN imgvr_hits_filt ON (tmp_table.cluster_id=imgvr_hits_filt.cluster_id) GROUP BY target_id, crispr_array) TO 'Target_to_repeat.tsv' (HEADER, DELIMITER '\t');
# And we remove the intermediary table
DROP TABLE tmp_table;