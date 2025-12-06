# Create a view, will be useful for later
CREATE OR REPLACE VIEW array_sample AS SELECT repeat_cluster, sra_run, ecosystem_sum, SUM(spacer_coverage) as total_coverage FROM array_tbl, sample_tbl, spacer_filt_tbl WHERE crispr_array==repeat_cluster AND spacer_filt_tbl.library==sample_tbl.library GROUP BY repeat_cluster, sra_run, ecosystem_sum;

# Making the counts for the paper 
## How many total spacers
SELECT COUNT(DISTINCT spacer_id) FROM spacer_filt_clusters;
## How many distinct spacers
SELECT COUNT(DISTINCT cluster_id) FROM spacer_filt_clusters;
## How many distinct samples
SELECT COUNT(DISTINCT sra_run) FROM array_sample;
## How many distinct repeats
SELECT COUNT(DISTINCT repeat_cluster) FROM array_sample;


# Create a file linking samples to a count of repeats - This only takes hq spacer, i.e. at least 5 distinct good-quality spacers
# All samples should have at least 1 repeat with high-quality spacers
COPY (SELECT * FROM (SELECT sra_run, COUNT(repeat_cluster) AS total_cluster FROM array_sample GROUP BY sra_run) AS tmp RIGHT JOIN sample_tbl ON (tmp.sra_run == sample_tbl.sra_run)) TO 'sample_to_array_counts.tsv' (HEADER, DELIMITER '\t', NULLSTR 'NA');

# Create a file linking repeats to count of samples
COPY (SELECT * FROM (SELECT repeat_cluster, COUNT(sra_run) AS total_sample FROM array_sample GROUP BY repeat_cluster) AS tmp RIGHT JOIN array_tbl ON (tmp.repeat_cluster == array_tbl.repeat_cluster)) TO 'array_to_sample_counts.tsv' (HEADER, DELIMITER '\t', NULLSTR 'NA');
# Same but this time counting the number of spacers for each repeat
COPY (SELECT lca_origin, SUM(total_spacers) as total_spacers FROM (SELECT crispr_array, count(*) as total_spacers FROM spacer_filt_tbl GROUP BY crispr_array) as tmp, array_tbl a WHERE a.repeat_cluster=tmp.crispr_array GROUP BY lca_origin) TO 'lca_origin_to_spacer_count.tsv' (HEADER, DELIMITER '\t');
# Also counting the number of repeats in each taxonomy category
COPY (SELECT lca_origin, count(*) AS total FROM array_tbl GROUP BY lca_origin) TO 'lca_origin_to_repeat_count-for-ref.tsv' (HEADER, DELIMITER '\t');
# And the same count, but including predicted CRISPR type information as well
COPY (SELECT type, lca_origin, count(*) AS total FROM array_tbl GROUP BY type, lca_origin) TO 'lca_origin_to_type.tsv' (HEADER, DELIMITER '\t');

# Look at distribution of individual spacers across arrays/repeats
COPY (SELECT n_array, COUNT(n_array) FROM (SELECT spacer_filt_clusters.cluster_id, COUNT(DISTINCT(tmp.crispr_array)) as n_array FROM (SELECT spacer_id, crispr_array, lca_origin, lca_genus FROM spacer_filt_tbl JOIN array_tbl ON array_tbl.repeat_cluster=spacer_filt_tbl.crispr_array) AS tmp JOIN spacer_filt_clusters ON spacer_filt_clusters.spacer_id=tmp.spacer_id GROUP BY spacer_filt_clusters.cluster_id) GROUP BY n_array) TO 'spacer_to_array_counts.tsv' (HEADER, DELIMITER '\t');

# Then look at only the one with a genus affiliation, and count how many distinct genera
COPY (SELECT n_array, n_genus, COUNT(*) FROM (SELECT spacer_filt_clusters.cluster_id, COUNT(DISTINCT(tmp.crispr_array)) as n_array, COUNT(DISTINCT(tmp.lca_genus)) as n_genus FROM (SELECT spacer_id, crispr_array, lca_origin, lca_genus FROM spacer_filt_tbl JOIN array_tbl ON array_tbl.repeat_cluster=spacer_filt_tbl.crispr_array WHERE lca_genus NOT LIKE '%g__unclassified' AND lca_genus NOT LIKE 'NA') AS tmp JOIN spacer_filt_clusters ON spacer_filt_clusters.spacer_id=tmp.spacer_id GROUP BY spacer_filt_clusters.cluster_id) GROUP BY n_genus, n_array) TO 'spacer_to_genus_counts.tsv' (HEADER, DELIMITER '\t');
