## In DuckDB
## Viruses of interest
CREATE TABLE tmp_taxa AS (SELECT * FROM read_csv('list_taxa.tsv',delim = '\t',header = true));
CREATE TABLE tmp_virus AS (SELECT * FROM (SELECT * FROM imgvr_info WHERE (host_from_taxo=='phage' OR host_from_taxo=='archaea')) AS tmp RIGHT JOIN tmp_taxa ON tmp.class==tmp_taxa.Relevant_taxa);
COPY (SELECT * FROM tmp_virus) TO 'selected_viruses.tsv' (HEADER, DELIMITER '\t');
## Spacers of interest
CREATE TABLE tmp_spacer AS (SELECT * FROM (SELECT uvig, cluster_id, hit_start, hit_end, hit_strand, n_mismatches FROM tmp_virus LEFT JOIN imgvr_hits_filt ON tmp_virus.uvig==imgvr_hits_filt.target_id) AS tmp LEFT JOIN spacer_filt_clusters ON tmp.cluster_id=spacer_filt_clusters.cluster_id);
SELECT * FROM tmp_spacer LIMIT 10;
COPY (SELECT * FROM tmp_spacer) TO 'selected_spacer_hits.tsv' (HEADER, DELIMITER '\t');
## Info on spacers of interest
CREATE TABLE list_spacer AS (SELECT DISTINCT(spacer_id) FROM tmp_spacer);
COPY (SELECT * FROM list_spacer LEFT JOIN spacer_filt_tbl ON list_spacer.spacer_id=spacer_filt_tbl.spacer_id LEFT JOIN array_tbl ON spacer_filt_tbl.crispr_array=array_tbl.repeat_cluster LEFT JOIN sample_tbl ON spacer_filt_tbl.library=sample_tbl.library) TO 'selected_spacer_info.tsv' (HEADER, DELIMITER '\t');
## Cleanup
DROP TABLE tmp_taxa;
DROP TABLE tmp_virus;
DROP TABLE tmp_spacer;
DROP TABLE list_spacer;