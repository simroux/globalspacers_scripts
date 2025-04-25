
## Intermediary table will be all spacers with at least 1 hit, list with their cluster_id, their library, and their array
CREATE TABLE spacer_hits_imgvr AS (SELECT spacer_filt_tbl.spacer_id as spacer_id, cluster_id, library, crispr_array, spacer_coverage FROM spacer_filt_tbl, (SELECT spacer_id, cluster_id FROM spacer_filt_clusters SEMI JOIN imgvr_hits_filt ON imgvr_hits_filt.cluster_id=spacer_filt_clusters.cluster_id) as tmp WHERE tmp.spacer_id=spacer_filt_tbl.spacer_id);

## How many non-redundant spacers have a hit:
## Look at a combined number for the abstract
select count(distinct(cluster_id)) from imgvr_hits_filt;
# 38,443,798 
select count(distinct(cluster_id)) from imgpr_hits_filt;
# 2,996,301
select count(distinct(cluster_id)) from imgpr_hits_filt semi join imgvr_hits_filt on imgpr_hits_filt.cluster_id=imgvr_hits_filt.cluster_id;
# 441,820
# So total: 38.4+3-0.440 -> 41 

## Export data for panel A
## Add info about hit / no hit based on imgvr_hits_filt
COPY (SELECT * FROM imgvr_info LEFT JOIN (SELECT DISTINCT(target_id), 1 as hit from imgvr_hits_filt) AS tmp ON imgvr_info.uvig==tmp.target_id) TO 'fig_hits_input_vr.tsv' (HEADER, DELIMITER '\t');
## Same with PR
COPY (SELECT * FROM imgpr_info LEFT JOIN (SELECT DISTINCT(target_id), 1 as hit from imgpr_hits_filt) AS tmp ON imgpr_info.full_plasmid_id==tmp.target_id) TO 'fig_hits_input_pr.tsv' (HEADER, DELIMITER '\t');

