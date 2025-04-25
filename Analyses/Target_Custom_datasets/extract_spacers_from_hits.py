#! /usr/bin/env python
import os
import sys
import argparse
import duckdb
import pandas as pd
import polars as pl


def fetch_arguments(parser):
	parser.add_argument('--in_file','-i', dest='list_hits', required=True, default='none',help='Table with the list of hits of interest (in a column "Spacer id")')
	parser.add_argument('--out_file','-o', dest='out_file', required=True, default='none',help='Output file')
	parser.add_argument('--path_db','-d', dest='path_db', required=True, default='none',help='Path to the duckdb file')

def get_sp_list(args):
	hit_df = pd.read_csv(args['list_hits'], sep='\t', low_memory=False)
	hit_df = hit_df.rename(columns={"Spacer id": "spid"})
	args["hit_df"] = hit_df

def get_spacer_info(args):
	with duckdb.connect(args["path_db"]) as con:
		df_tmp = args["hit_df"]
		req = f"SELECT * FROM (SELECT * FROM (SELECT * FROM df_tmp, spacer_filt_clusters WHERE spid=cluster_id) AS tmp LEFT JOIN spacer_tbl ON tmp.spacer_id=spacer_tbl.spacer_id) AS tmp_2 LEFT JOIN array_tbl ON tmp_2.crispr_array=array_tbl.repeat_cluster"
		print(f"{req}")
		# spacer_df = con.sql(req).df()
		spacer_df = con.sql(req).pl()
		print(f"{spacer_df}")
		args["spacer_df"] = spacer_df

def main():
	parser = argparse.ArgumentParser()
	fetch_arguments(parser)
	args = vars(parser.parse_args())
	print(f"{args['out_file']}")
	print(f".. load the list of spacer ids from {args['list_hits']}")
	get_sp_list(args)
	print(f".. get full complement of spacers from these samples")
	get_spacer_info(args)
	## Export
	args["spacer_df"].write_csv(args["out_file"], separator='\t')

if __name__ == "__main__":
	output = main()

