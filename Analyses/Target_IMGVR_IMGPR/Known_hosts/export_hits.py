#! /usr/bin/env python
import os
import sys
import argparse
import duckdb
import polars as pl

def fetch_arguments(parser):
    parser.add_argument('--uvig_list','-u', dest='in_file', required=True, default='none',help='file with the list of phages with a known host')
    parser.add_argument('--duckdb','-d', dest='path_db', required=True, default='none',help='path to the duckdb database')

def load_uvig_list(args):
    df_uvig = pl.read_csv(args["in_file"], separator='\t', low_memory=False)
    args["df_uvig"] = df_uvig

def get_hits(args):
    with duckdb.connect(args["path_db"],read_only = True) as con:
        df_tmp = args["df_uvig"]
        req = f"SELECT df_tmp.*, imgvr_hits_filt.* FROM imgvr_hits_filt, df_tmp WHERE df_tmp.cluster_id=imgvr_hits_filt.cluster_id"
        req = f"select target_id, n_hits, array_tbl.* from (select target_id, crispr_array, count(distinct(tmp.cluster_id)) as n_hits FROM (select cluster_id, target_id from imgvr_hits_filt, df_tmp where imgvr_hits_filt.target_id=df_tmp.uvig_id) AS tmp JOIN spacer_hits_imgvr ON tmp.cluster_id=spacer_hits_imgvr.cluster_id GROUP BY (target_id, crispr_array)) AS tmp2 JOIN array_tbl ON tmp2.crispr_array=array_tbl.repeat_cluster;"
        # req =
        args["df_hits"] = con.sql(req).pl()

def main():
    parser = argparse.ArgumentParser()
    fetch_arguments(parser)
    args = vars(parser.parse_args())
    print(f".. load the list of relevant uvigs from {args['in_file']}")
    load_uvig_list(args)
    print(f".. get full list of hits with repeat id for these uvigs")
    args["out_file_hits"] = "Phages_with_known_hosts-all_hit_counts_to_repeats.tsv"
    get_hits(args)
    args["df_hits"].write_csv(args["out_file_hits"], separator='\t')

if __name__ == "__main__":
    output = main()


