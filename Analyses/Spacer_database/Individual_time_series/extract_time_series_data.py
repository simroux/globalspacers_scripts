#! /usr/bin/env python
import os
import sys
import argparse
import duckdb
# import pandas as pd
import polars as pl

def fetch_arguments(parser):
    parser.add_argument('--series_id','-id', dest='series_id', required=True, default='none',help='ID of the series of interest (e.g. PRJEB62690_025)')

def load_lib_list(args):
    df_sample = pl.read_csv(args["in_file"], separator='\t', low_memory=False)
    args["df_sample"] = df_sample[["bioproject","subject_id","date","library"]]

def get_spacer_table(args):
    with duckdb.connect(args["path_db"],read_only = True) as con:
        df_tmp = args["df_sample"]
        req = f"SELECT tmp.*, cluster_id FROM (SELECT spacer_filt_tbl.*, bioproject, subject_id, date FROM spacer_filt_tbl, df_tmp WHERE df_tmp.library=spacer_filt_tbl.library) AS tmp LEFT JOIN spacer_clusters ON tmp.spacer_id=spacer_clusters.spacer_id"
        print(f"{req}")
        spacer_df = con.sql(req).pl()
        print(f"{spacer_df}")
        args["spacer_df"] = spacer_df

def get_hits(args):
    with duckdb.connect(args["path_db"],read_only = True) as con:
        df_tmp = args["spacer_df"]
        req = f"SELECT df_tmp.*, imgvr_hits_filt.* FROM imgvr_hits_filt, df_tmp WHERE df_tmp.cluster_id=imgvr_hits_filt.cluster_id"
        print(f"{req}")
        df_hits_one = con.sql(req).pl()
        req = f"SELECT df_tmp.*, imgpr_hits_filt.* FROM imgpr_hits_filt, df_tmp WHERE df_tmp.cluster_id=imgpr_hits_filt.cluster_id"
        print(f"{req}")
        df_hits_two = con.sql(req).pl()
        req = f"SELECT df_tmp.*, imgvr_hits_extra_filt.* FROM imgvr_hits_extra_filt, df_tmp WHERE df_tmp.cluster_id=imgvr_hits_extra_filt.cluster_id"
        print(f"{req}")
        df_hits_three = con.sql(req).pl()
        req = f"SELECT df_tmp.*, imgpr_hits_extra_filt.* FROM imgpr_hits_extra_filt, df_tmp WHERE df_tmp.cluster_id=imgpr_hits_extra_filt.cluster_id"
        print(f"{req}")
        df_hits_four = con.sql(req).pl()
        args["df_hits"] = pl.concat([df_hits_one,df_hits_two,df_hits_three,df_hits_four])

def main():
    parser = argparse.ArgumentParser()
    fetch_arguments(parser)
    args = vars(parser.parse_args())
    args["data_dir"] = "Data/Spacer_db/"
    args["path_db"] = "Data/global_crispr_db.duckdb"
    args["wdir"] = "./"
    ## 
    args["in_file"] = os.path.join(args["wdir"],args["series_id"],args["series_id"]+"_samples.tsv")
    print(f".. load the list of relevant samples from {args['in_file']} corresponding to {args['series_id']}")
    load_lib_list(args)
    ## 
    print(f".. get full complement of spacers from these libraries")
    args["out_file_spacers"] = os.path.join(args["wdir"],args["series_id"],args["series_id"]+"_filt_spacers.tsv")
    get_spacer_table(args)
    args["spacer_df"].write_csv(args["out_file_spacers"], separator='\t')
    ## 
    print(f" .. get all uvig and plasmid hits")
    args["out_file_hits"] = os.path.join(args["wdir"],args["series_id"],args["series_id"]+"_filt_hits.tsv")
    get_hits(args)
    args["df_hits"].write_csv(args["out_file_hits"], separator='\t')

if __name__ == "__main__":
    output = main()

