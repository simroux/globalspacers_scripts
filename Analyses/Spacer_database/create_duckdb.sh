################## CREATE DATABASE ####################
### Note: you may not be able to run all these commands in the same session, and may need to flush the memory from time to time by quitting DuckDb, then reloading the db
duckdb global_crispr_db.duckdb

### Load samples
CREATE TABLE sample_tbl AS FROM read_csv('Data/Spacer_db/Runs_to_ecosystem_and_sequencing_and_study_for_db-Jul28-24.tsv', header = true, delim='\t');
CREATE INDEX lib_id_ix ON sample_tbl (library);
DESCRIBE sample_tbl;

### Load array info
CREATE TABLE array_tbl AS FROM read_csv('Data/Spacer_db/Array_info_filtered_for_db-Oct24-25.tsv', header = true, delim='\t');
CREATE UNIQUE INDEX array_id_ix ON array_tbl (repeat_cluster);
DESCRIBE array_tbl;

### Load spacers
CREATE TABLE spacer_tbl AS FROM read_csv('Data/Spacer_db/All_spacers_info_filtered-Jul19-24.tsv', delim = '\t', header = true);
ALTER TABLE spacer_tbl ALTER hq TYPE TINYINT;
CREATE INDEX spacer_id_idx ON spacer_tbl (spacer_id); 
## Note: this is not unique, because some spacer ids are duplicates (when found in multiple arrays)
DESCRIBE spacer_tbl;

### Load spacer clusters
CREATE TABLE spacer_clusters AS FROM read_csv('Data/Spacer_db/All_spacers_info_filtered_clusters-Jul19-24.tsv', delim = '\t', header = true);
CREATE INDEX cls_spacer_id_idx ON spacer_clusters (spacer_id);
CREATE INDEX cls_cluster_id_idx ON spacer_clusters (cluster_id);
DESCRIBE spacer_clusters;

### Create additional tables with only hq spacers, so that we don't have to redo this (massive) join every time
# Note -> Not ideal from a space perspective, but super useful from a memory/time perspective later
CREATE TABLE spacer_filt_tbl AS (SELECT * FROM spacer_tbl WHERE hq==1);
CREATE TABLE spacer_filt_clusters AS (SELECT cluster_id, spacer_id FROM spacer_clusters SEMI JOIN spacer_filt_tbl ON spacer_clusters.spacer_id = spacer_filt_tbl.spacer_id);

