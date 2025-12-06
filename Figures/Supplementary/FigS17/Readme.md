# How-to - reconstruct panels of Fig. S17
First, we need to export some information from the DuckDB spacer database:
```
COPY (select DISTINCT crispr_array, lca_genus, s.library, ecosystem_sum FROM (select t.crispr_array, a.lca_genus, t.library from spacer_filt_tbl t, array_tbl a WHERE lca_genus NOT LIKE 'NA' AND a.repeat_cluster=t.crispr_array) as tmp, sample_tbl s WHERE tmp.library=s.library) TO 'taxon_to_ecosystem.tsv' (HEADER, DELIMITER '\t');
```

Then, we use the script `./get_repeat_cotargeting.pl run` to generate the file `cotargeting_vs_codetection_min100.tsv`. This file will then be the input in `plot_codetection.R` to prepare the plot.

