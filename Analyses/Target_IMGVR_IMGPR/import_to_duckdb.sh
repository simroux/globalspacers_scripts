## In DuckDB:
## Import all hits to IMG/VR
## and create a filtered table as well here
CREATE TABLE imgvr_hits AS FROM read_csv('../../Data/Results_nr_imgvrhq_all_best_matches-Sept13-24.tsv', delim = '\t', header = true);
CREATE TABLE imgvr_hits_filt AS (SELECT * FROM imgvr_hits SEMI JOIN spacer_filt_clusters ON imgvr_hits.cluster_id = spacer_filt_clusters.cluster_id);
## We also import the IMG/VR information on each UVIG, and we checked we had the id that match (UVig_id|taxon|scaffold)
CREATE TABLE imgvr_info AS FROM read_csv('../../Data/IMGVR_sequence_information_Oct17.tsv', delim = '\t', header = true);
# Check if we have any hit targets with no information
SELECT target_id FROM imgvr_hits_filt ANTI JOIN imgvr_info ON imgvr_hits_filt.target_id = imgvr_info.uvig;
# Confirmed, 0

###### Import all the hits to IMG/PR
CREATE TABLE imgpr_hits AS FROM read_csv('../../Data/Results_nr_imgprhq_all_best_matches-Sept13-24.tsv', delim = '\t', header = true);
CREATE TABLE imgpr_hits_filt AS (SELECT * FROM imgpr_hits SEMI JOIN spacer_filt_clusters ON imgpr_hits.cluster_id = spacer_filt_clusters.cluster_id);
## We also import the IMG/VR information on each UVIG, and we checked we had the id that match (UVig_id|taxon|scaffold)
CREATE TABLE imgpr_info AS FROM read_csv('../../Data/IMGPR_sequence_information_Aug26.tsv', delim = '\t', header = true);
# Check if we have any hit targets with no information
SELECT target_id FROM imgpr_hits_filt ANTI JOIN imgpr_info ON imgpr_hits_filt.target_id = imgpr_info.full_plasmid_id;
## Confirmed, 0
