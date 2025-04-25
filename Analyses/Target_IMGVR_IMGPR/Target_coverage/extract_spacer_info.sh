## For coverage analysis, export all hits to high-quality viruses
# Create a table of hits restricted to only HQ UViGs (for coverage analysis)
CREATE TABLE tmp_table_1 AS (SELECT cluster_id, target_id, votu, hit_start, hit_end, hit_strand, n_mismatches FROM (SELECT cluster_id, target_id, hit_start, hit_end, hit_strand, n_mismatches FROM imgvr_hits_filt) AS tmp1 JOIN (SELECT uvig, votu FROM imgvr_info WHERE quality='High-quality' OR quality='Reference') AS tmp2 ON tmp1.target_id=tmp2.uvig);
COPY tmp_table_1 TO 'All_hits_vr_hq.tsv' (HEADER, DELIMITER '\t');
# Get the list of spacers associated with these hits and get all the info
CREATE TABLE tmp_table_2 AS (SELECT * FROM spacer_clusters SEMI JOIN tmp_table_1 USING (cluster_id));
COPY (SELECT * FROM tmp_table_2 JOIN spacer_filt_tbl USING (spacer_id)) TO 'All_spacers_vr_hq.tsv' (HEADER, DELIMITER '\t');
# 
DROP TABLE tmp_table_1;
DROP TABLE tmp_table_2;
## Proceed with the same approach for plasmids
# Create a table of hits restricted to only near-complete plasmids
CREATE TABLE tmp_table_1 AS (SELECT cluster_id, target_id, ptu, hit_start, hit_end, hit_strand, n_mismatches FROM (SELECT cluster_id, target_id, hit_start, hit_end, hit_strand, n_mismatches FROM imgpr_hits_filt) AS tmp1 JOIN (SELECT full_plasmid_id, ptu FROM imgpr_info WHERE putatively_complete='Yes') AS tmp2 ON tmp1.target_id=tmp2.full_plasmid_id);
COPY tmp_table_1 TO 'All_hits_pr_complete.tsv' (HEADER, DELIMITER '\t');
# Get the list of spacers associated with these hits and get all the info
CREATE TABLE tmp_table_2 AS (SELECT * FROM spacer_clusters SEMI JOIN tmp_table_1 USING (cluster_id));
COPY (SELECT * FROM tmp_table_2 JOIN spacer_filt_tbl USING (spacer_id)) TO 'All_spacers_pr_complete.tsv' (HEADER, DELIMITER '\t');
# 
DROP TABLE tmp_table_1;
DROP TABLE tmp_table_2;

## For analysis of viruses targeted by multiple taxa, we export similar information for viruses not already processed as high-quality
## The list can be generated using scripts provided as part of Figure 5 preparation, but we also provide the list directly in Additional_data
CREATE OR REPLACE TABLE multitaxa_uvig_list AS FROM read_csv('../../../Data/Additional_data/List_multiclass_uvigs_nohq.txt', delim = '\t', header = true);
# Create a table of hits restricted to these viruses
CREATE TABLE tmp_table_1 AS (SELECT cluster_id, target_id, votu, hit_start, hit_end, hit_strand, n_mismatches FROM (SELECT cluster_id, target_id, hit_start, hit_end, hit_strand, n_mismatches FROM imgvr_hits_filt, multitaxa_uvig_list WHERE multitaxa_uvig_list.uvig=imgvr_hits_filt.target_id) AS tmp JOIN imgvr_info ON imgvr_info.uvig=tmp.target_id);
COPY tmp_table_1 TO 'All_hits_vr_additional_multitaxa.tsv' (HEADER, DELIMITER '\t');
# Get the list of spacers associated with these hits and get all the info
CREATE TABLE tmp_table_2 AS (SELECT * FROM spacer_clusters SEMI JOIN tmp_table_1 USING (cluster_id));
COPY (SELECT * FROM tmp_table_2 JOIN spacer_filt_tbl USING (spacer_id)) TO 'All_spacers_vr_additional_multitaxa.tsv' (HEADER, DELIMITER '\t');
# Remove tmp tables
DROP TABLE multitaxa_uvig_list;
DROP TABLE tmp_table_1;
DROP TABLE tmp_table_2;