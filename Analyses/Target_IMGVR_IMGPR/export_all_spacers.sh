## In DuckDB:
COPY (SELECT DISTINCT cluster_id, spacer_sequence FROM spacer_tbl sp, spacer_clusters sc WHERE sc.spacer_id==sp.spacer_id AND sp.hq==1) TO 'nr_spacers_hq.tsv' (HEADER, DELIMITER '\t');
## Then in bash / the terminal:
cat nr_spacers_hq.tsv | grep -v cluster_id | gawk -F '\t' '{print ">"$0}' | tr '\t' '\n' > nr_spacers_hq.fna
