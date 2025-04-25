#! /usr/bin/env python
import os
import sys
import argparse
import duckdb
import pandas as pd
import polars as pl

def fetch_arguments(parser):
	parser.add_argument('--uvig_id','-id', dest='uvig_id', required=True, default='none',help='ID of the UViG of interest (e.g. IMGVR_UViG_3300007362_000033)')


def guess_full_uvig(args):
	imgvr_info_file = "../../../Data/Additional_data/translate_uvig.tsv"
	req = f"SELECT uvig_id_full FROM '{imgvr_info_file}' WHERE uvig_id=='{args['uvig_id']}'"
	print(f"{req}")
	toto = duckdb.sql(req).df()
	args['uvig_id_full'] = toto["uvig_id_full"].values[0]
	print(f"{args['uvig_id']} => {args['uvig_id_full']}")

def get_spacer_hits(args):
	with duckdb.connect(args["path_db"]) as con:
		req = f"SELECT * FROM imgvr_hits_filt, spacer_filt_clusters, spacer_filt_tbl, sample_tbl, array_tbl WHERE target_id='{args['uvig_id_full']}' AND imgvr_hits_filt.cluster_id=spacer_filt_clusters.cluster_id AND spacer_filt_clusters.spacer_id=spacer_filt_tbl.spacer_id AND spacer_filt_tbl.library=sample_tbl.library AND array_tbl.repeat_cluster=spacer_filt_tbl.crispr_array"
		print(f"{req}")
		# spacer_df = con.sql(req).df()
		spacer_df = con.sql(req).pl()
		print(f"{spacer_df}")
		args["spacer_df"] = spacer_df

def add_nvotu_cluster(args):
	## Polars
	df_nvotu = pl.read_csv(args["nvotu_cluster"], separator='\t', low_memory=False)
	df1 = args["spacer_df"]
	df3 = df1.join(df_nvotu, on='cluster_id', how='left')
	# df3["n_votu"].fill_null(1)
	args["spacer_df"] = df3

def add_array_size(args):
	df_arraysize = pl.read_csv(args["alphadiv"], separator='\t', low_memory=False)
	df1 = args["spacer_df"]
	df3 = df1.join(df_arraysize[["array","sample","n_spacers","max_cover"]], left_on=['crispr_array','sra_run'], right_on=['array','sample'], how='left')
	args["spacer_df"] = df3
	
def add_common_spacers(args):
	### Polars
	df_common = pl.read_csv(args["common"], separator='\t', low_memory=False)
	# sp_list = ["Sp_cl_000487759904","Sp_cl_000763259640","Sp_cl_000547251874","Sp_cl_000000000056","Sp_cl_000000000053"]
	sp_list = args["spacer_df"]["cluster_id"].unique()
	df_common = df_common.filter(pl.col('spacer').is_in(sp_list))
	## select columns of interest
	df_common = df_common.select(["spacer","list_common_sets"])
	## transform the list_common_sets in proper lists 
	df_common = df_common.with_columns(pl.col("list_common_sets").str.split("|").alias("list_common_sets"))
	## explode
	df_common = df_common.explode("list_common_sets")
	## get the array and sample separately
	df_common = df_common.with_columns(pl.col("list_common_sets").str.split_exact(",", 1).struct.rename_fields(["array","sample"]).alias("fields")).unnest("fields")
	## remove column we don't need anymore, and add one with yes
	df_common = df_common.select(["spacer","array","sample"])
	df_common = df_common.with_columns(pl.lit("yes").alias("is_common"))
	print(f"{df_common}")
	#### Merging with main data frame
	df1 = args["spacer_df"]
	df3 = df1.join(df_common[["spacer","array","sample","is_common"]], left_on=['cluster_id','crispr_array','sra_run'], right_on=['spacer','array','sample'], how='left')
	df3 = df3.with_columns(pl.col("is_common").fill_null("no"))
	args["spacer_df"] = df3

def mergeOverlap(arr):
    # Sort intervals based on start values
    arr.sort()
    res = []
    res.append(arr[0])
    for i in range(1, len(arr)):
        last = res[-1]
        curr = arr[i]
        # If current interval overlaps with the last merged
        # interval, merge them 
        if curr[0] <= last[1]:
            last[1] = max(last[1], curr[1])
        else:
            res.append(curr)
    return res

def get_cover(vec):
    ## merge all overlapping intervals
    list_start_end = vec.with_columns(concat_list=pl.concat_list("hit_start", "hit_end"))['concat_list'].to_list()
    res = mergeOverlap(list_start_end)
    ## get length
    total_cover = 0
    for interval in res:
        total_cover = total_cover + (interval[1]-interval[0]+1)
        # print(interval[0], interval[1])
    array = vec.select(pl.first("crispr_array"))["crispr_array"].to_list()
    df = pl.DataFrame({"crispr_array": array,"all_samples_cover": [total_cover],})
    return df

def get_cover_bysample(vec):
    ## merge all overlapping intervals
    list_start_end = vec.with_columns(concat_list=pl.concat_list("hit_start", "hit_end"))['concat_list'].to_list()
    res = mergeOverlap(list_start_end)
    ## get length
    total_cover = 0
    for interval in res:
        total_cover = total_cover + (interval[1]-interval[0]+1)
    array = vec.select(pl.first("crispr_array"))["crispr_array"].to_list()
    sra_run = vec.select(pl.first("sra_run"))["sra_run"].to_list()
    df = pl.DataFrame({"crispr_array": array,"sra_run": sra_run,"individual_sample_cover": [total_cover],})
    return df

def add_coverage_all_samples(args):
    df_cover_all = args["spacer_df"].group_by("crispr_array").map_groups(lambda x: get_cover(x))
    # print(f"{toto}")
    df3 = args["spacer_df"].join(df_cover_all[["crispr_array","all_samples_cover"]], left_on=['crispr_array'], right_on=['crispr_array'], how='left')
    args["spacer_df"] = df3

def add_coverage_by_sample(args):
    df_cover_all = args["spacer_df"].group_by("crispr_array","sra_run").map_groups(lambda x: get_cover_bysample(x))
    # print(f"{toto}")
    df3 = args["spacer_df"].join(df_cover_all[["crispr_array","sra_run","individual_sample_cover"]], left_on=['crispr_array','sra_run'], right_on=['crispr_array','sra_run'], how='left')
    args["spacer_df"] = df3


def main():
	parser = argparse.ArgumentParser()
	fetch_arguments(parser)
	args = vars(parser.parse_args())
	args["out_file"] = os.path.join(f"{args['uvig_id']}_spacer_hits_full.tsv")
	# args["out_file_clean"] = os.path.join(f"{args['uvig_id']}_spacer_hits_clean.tsv")
	args["data_dir"] = "../../../Data/Spacer_db/"
	args["path_db"] = "../../../Analyses/Spacer_database/global_crispr_db.duckdb"
	args["nvotu_cluster"] = "/../../../Analyses/Target_IMGVR_IMGPR/n_votu_by_cluster.tsv"
	args["alphadiv"] = "../../../Analyses/Spacer_database/spacer_sets_alphadiv.tsv"
	args["common"] = "../../../Analyses/Spacer_database/Common_spacer_list.tsv"
	print(f"{args['out_file']}")
	## 
	print(f".. load the full uvig_id corresponding to {args['uvig_id']}")
	guess_full_uvig(args)
	## 
	print(f".. get full complement of spacers hitting this uvig")
	get_spacer_hits(args)
	## 
	print(f".. calculate coverage per array across all samples")
	add_coverage_all_samples(args)
	## 
	print(f".. calculate coverage per array for each sample")
	add_coverage_by_sample(args)
    ## 
	print(f".. load number of votu by spacer cluster")
	add_nvotu_cluster(args)
	## 
	print(f".. load information about array size")
	add_array_size(args)
	##
	print(f".. load information about common spacers")
	add_common_spacers(args)
	## Export
	args["spacer_df"].write_csv(args["out_file"], separator='\t')
      
if __name__ == "__main__":
	output = main()

