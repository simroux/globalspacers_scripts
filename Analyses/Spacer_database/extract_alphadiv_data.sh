CREATE OR REPLACE VIEW spacer_cover AS SELECT cluster_id, SUM(spacer_coverage) as total_coverage, crispr_array, sra_run FROM spacer_filt_clusters sc, spacer_filt_tbl st, sample_tbl as sp WHERE st.spacer_id==sc.spacer_id AND st.library==sp.library GROUP BY cluster_id,crispr_array,sra_run;
## Simple export first, listing the distribution of each (non-redundant) spacer
COPY (SELECT cluster_id, crispr_array, count(sra_run) as n_sample, sum(singleton) as n_singleton, list(sra_run) as list_sample FROM (SELECT cluster_id, crispr_array, sra_run, CASE WHEN total_coverage > 1 THEN 0 ELSE 1 END as singleton FROM spacer_cover) GROUP BY cluster_id, crispr_array ORDER BY crispr_array) TO 'spacer_clusters_and_metadata.tsv' (HEADER, DELIMITER '\t');
## Next make an export a bit more complicated that acutally includes the coverage of each spacer in each sample
## We do this in two steps to keep the query "not too overwhelming" in terms of memory
## For sets with more than 10 clusters (i.e. more than 10 distinct spacers)
CREATE OR REPLACE VIEW selected_sets_for_alphadiv AS SELECT * FROM (SELECT crispr_array, sra_run, count(distinct(cluster_id)) as total_cluster FROM spacer_cover GROUP BY crispr_array, sra_run) WHERE total_cluster>=10;
COPY(select sc.cluster_id, sc.total_coverage, sc.crispr_array, sc.sra_run from spacer_cover sc right join selected_sets_for_alphadiv sa on (sa.crispr_array=sc.crispr_array AND sa.sra_run=sc.sra_run) ORDER BY sc.crispr_array, sc.sra_run) TO 'spacer_clusters_and_metadata_for_alphadiv.tsv' (HEADER, DELIMITER '\t');
## For the ones below 10 clusters
CREATE OR REPLACE VIEW selected_sets_for_alphadiv AS SELECT * FROM (SELECT crispr_array, sra_run, count(distinct(cluster_id)) as total_cluster FROM spacer_cover GROUP BY crispr_array, sra_run) WHERE total_cluster<10;
COPY(select sc.cluster_id, sc.total_coverage, sc.crispr_array, sc.sra_run from spacer_cover sc right join selected_sets_for_alphadiv sa on (sa.crispr_array=sc.crispr_array AND sa.sra_run=sc.sra_run) ORDER BY sc.crispr_array, sc.sra_run) TO 'spacer_clusters_and_metadata_for_alphadiv_md10.tsv' (HEADER, DELIMITER '\t');
## Also get a count of ecosystems per spacer cluster, which will be used when looking at spacer with/without a hit in IMG/VR / PR (and when comparing the ecosystem of the spacer and of the target)
COPY (select cluster_id, ecosystem_sum, count(*) as total FROM (select cluster_id, library FROM spacer_filt_clusters JOIN spacer_filt_tbl ON spacer_filt_clusters.spacer_id=spacer_filt_tbl.spacer_id) AS tmp_1 JOIN sample_tbl ON tmp_1.library=sample_tbl.library GROUP BY cluster_id, ecosystem_sum) TO 'spacer_clusters_to_ecosystem.tsv' (HEADER, DELIMITER '\t');

