## In duckdb
## Import all extra hits to IMG/VR
CREATE TABLE imgvr_hits_extra AS FROM read_csv('../../../Data/Spacer_db/Additional_hits_imgvr-Dec13.tsv', delim = '\t', header = true);
CREATE TABLE imgvr_hits_extra_filt AS (SELECT * FROM imgvr_hits_extra SEMI JOIN spacer_filt_clusters ON imgvr_hits_extra.cluster_id = spacer_filt_clusters.cluster_id);
## For coverage analysis, export all hits to high-quality viruses
CREATE TABLE tmp_table_1 AS (SELECT cluster_id, target_id, votu, hit_start, hit_end, hit_strand, n_mismatches FROM (SELECT cluster_id, target_id, hit_start, hit_end, hit_strand, n_mismatches FROM imgvr_hits_extra_filt) AS tmp1 JOIN (SELECT uvig, votu FROM imgvr_info WHERE quality='High-quality' OR quality='Reference') AS tmp2 ON tmp1.target_id=tmp2.uvig); ### Note - we should already have only the high-quality, but we still need the votu 
COPY tmp_table_1 TO 'All_hits_vr_hq_2-3mis.tsv' (HEADER, DELIMITER '\t');
# Get the list of spacers associated with these hits and get all the info
CREATE TABLE tmp_table_2 AS (SELECT * FROM spacer_clusters SEMI JOIN tmp_table_1 USING (cluster_id));
COPY (SELECT * FROM tmp_table_2 JOIN spacer_filt_tbl USING (spacer_id)) TO 'All_spacers_vr_hq_2-3mis.tsv' (HEADER, DELIMITER '\t');
# 
DROP TABLE tmp_table_1;
DROP TABLE tmp_table_2;

## Same for plasmids
## Import all extra hits to IMG/PR
CREATE TABLE imgpr_hits_extra AS FROM read_csv('../../Data/Spacer_db/Additional_hits_imgpr-Dec13.tsv', delim = '\t', header = true);
CREATE TABLE imgpr_hits_extra_filt AS (SELECT * FROM imgpr_hits_extra SEMI JOIN spacer_filt_clusters ON imgpr_hits_extra.cluster_id = spacer_filt_clusters.cluster_id);
## For coverage analysis, export all hits to high-quality viruses
CREATE TABLE tmp_table_1 AS (SELECT cluster_id, target_id, ptu, hit_start, hit_end, hit_strand, n_mismatches FROM (SELECT cluster_id, target_id, hit_start, hit_end, hit_strand, n_mismatches FROM imgpr_hits_extra_filt) AS tmp1 JOIN (SELECT full_plasmid_id, ptu FROM imgpr_info WHERE putatively_complete='Yes') AS tmp2 ON tmp1.target_id=tmp2.full_plasmid_id); 
COPY tmp_table_1 TO 'All_hits_pr_hq_2-3mis.tsv' (HEADER, DELIMITER '\t');
# Get the list of spacers associated with these hits and get all the info
CREATE TABLE tmp_table_2 AS (SELECT * FROM spacer_clusters SEMI JOIN tmp_table_1 USING (cluster_id));
COPY (SELECT * FROM tmp_table_2 JOIN spacer_filt_tbl USING (spacer_id)) TO 'All_spacers_pr_hq_2-3mis.tsv' (HEADER, DELIMITER '\t');
# 
DROP TABLE tmp_table_1;
DROP TABLE tmp_table_2;

## Same for multitaxa-targeted viruses
###### First: export the 0 and 1 mismatch hits if not already exported for coverage info
 ## If not already loaded:
CREATE OR REPLACE TABLE multitaxa_uvig_list AS FROM read_csv('../../../Data/Additional_data/List_multiclass_uvigs_nohq.txt', delim = '\t', header = true);
###### Import all extra hits to IMG/VR
CREATE TABLE imgvr_hits_multitaxa_tmp AS FROM read_csv('../../Data/Spacer_db/Additional_hits_imgvr-multitaxa-Dec17.tsv', delim = '\t', header = true);
CREATE TABLE imgvr_hits_multitaxa_filt AS (SELECT * FROM imgvr_hits_multitaxa_tmp SEMI JOIN spacer_filt_clusters ON imgvr_hits_multitaxa_tmp.cluster_id = spacer_filt_clusters.cluster_id);
DROP TABLE imgvr_hits_multitaxa_tmp;
## For coverage analysis, export all hits to relevant viruses
CREATE TABLE tmp_table_1 AS (SELECT cluster_id, target_id, votu, hit_start, hit_end, hit_strand, n_mismatches FROM (SELECT cluster_id, target_id, hit_start, hit_end, hit_strand, n_mismatches FROM imgvr_hits_multitaxa_filt) AS tmp1 JOIN (SELECT uvig, votu FROM imgvr_info) AS tmp2 ON tmp1.target_id=tmp2.uvig); ### Note - we should already have only the high-quality, but we still need the votu 
COPY tmp_table_1 TO 'All_hits_vr_additional_multitaxa_2-3mis.tsv' (HEADER, DELIMITER '\t');
# Get the list of spacers associated with these hits and get all the info
CREATE TABLE tmp_table_2 AS (SELECT * FROM spacer_clusters SEMI JOIN tmp_table_1 USING (cluster_id));
COPY (SELECT * FROM tmp_table_2 JOIN spacer_filt_tbl USING (spacer_id)) TO 'All_spacers_vr_additional_multitaxa_2-3mis.tsv' (HEADER, DELIMITER '\t');
# 
DROP TABLE tmp_table_1;
DROP TABLE tmp_table_2;
