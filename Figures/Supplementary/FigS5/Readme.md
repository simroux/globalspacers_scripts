# Input file generation
The input files for this supplementary figure are "genus_to_ecosystem_counts.tsv" and "type_vs_ecosystem.txt". These can be generated from the duckdb spacer database as follows:
```
CREATE OR REPLACE VIEW array_sample AS SELECT repeat_cluster, sra_run, ecosystem_sum, SUM(spacer_coverage) as total_coverage FROM array_tbl, sample_tbl, spacer_filt_tbl WHERE crispr_array==repeat_cluster AND spacer_filt_tbl.library==sample_tbl.library GROUP BY repeat_cluster, sra_run, ecosystem_sum;

COPY (SELECT ecosystem_sum, lca_genus, count(tmp.repeat_cluster) as total_clusters FROM (SELECT repeat_cluster, ecosystem_sum from array_sample) AS tmp JOIN array_tbl ON (tmp.repeat_cluster == array_tbl.repeat_cluster) WHERE lca_genus <> 'NA' GROUP BY (ecosystem_sum, lca_genus)) TO 'genus_to_ecosystem_counts.tsv' (HEADER, DELIMITER '\t');

COPY (SELECT type, ecosystem_sum, count(*) as n_occurrences from array_sample group by type, ecosystem_sum) TO 'type_vs_ecosystem.txt' (HEADER, DELIMITER '\t', NULLSTR 'NA');

```  


# Use R to generate the figure
`plot_fig_ecosystem_vs_repeat_taxo.R`
To generate the figures, in R.
