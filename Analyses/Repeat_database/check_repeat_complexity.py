#!/usr/bin/env python
import argparse
import csv
import os
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.lcc import lcc_simp
import kcounter

def main():
	parser = argparse.ArgumentParser()
	# parser.add_argument("-d", "--repeat_db_dir", dest='db_dir', required=True, default='none', help="path to the repeat database folder")
	# parser.add_argument("-w", "--wdir", dest='wdir', required=True, default='none', help="output directory (will be created if it does not exist)")
	# parser.add_argument("-t", "--n_threads", dest='n_threads', required=False, default=2, help="number of threads to use (default 2)")
	args = vars(parser.parse_args())
	set_default(args)
	# 
	if not (os.path.exists(args["out_file"])): ## Creating output file if needed
		with open(args["out_file"], "w") as output_handle_tab, open(args["out_file_detail"], "w") as output_handle_tab_bis:
			out_tab_writer = csv.writer(output_handle_tab, delimiter='\t',quotechar='"', quoting=csv.QUOTE_MINIMAL, lineterminator="\n")
			out_tab_writer.writerow(["Repeat","Lcc","Kmer","Final_flag","Sequence"])
			out_detail_tab_writer = csv.writer(output_handle_tab_bis, delimiter='\t',quotechar='"', quoting=csv.QUOTE_MINIMAL, lineterminator="\n")
			out_detail_tab_writer.writerow(["Seq","Kmer","Max","Sum","Ratio","Flag_low_complexity","Counts"])
			# seq,i,max(counts),sum(counts),ratio,counts
			for seq_record in SeqIO.parse(args["input_repeats"], "fasta"):
				flag_lcc = test_lcc(seq_record.seq)
				flag_kmer = test_kmer(str(seq_record.seq),out_detail_tab_writer)
				final_flag = ""
				if flag_lcc != "no" or flag_kmer != "no":
					final_flag = "Low_complexity"
				print(f"Result\t{seq_record.id}\t{flag_lcc}\t{flag_kmer}\t{seq_record.seq}")
				out_tab_writer.writerow([seq_record.id,flag_lcc,flag_kmer,final_flag,seq_record.seq])
	df_complex = pd.read_csv(args["out_file"],sep="\t")
	print(f"{df_complex}")
	df_info = pd.read_csv(args["ref_file"],sep="\t")
	print(f"{df_info}")
	df_info = df_info.merge(df_complex, left_on="Cluster", right_on="Repeat", how="left")
	df_info.to_csv("Array_info_with_complexity.tsv",index=False, sep='\t', na_rep='NA')

	

def test_lcc(seq):
	seq_lcc = lcc_simp(seq)
	flag = "no"
	if seq_lcc < 0.8:
		flag = "low_complexity_lcc"
	return flag

def test_kmer(seq,out):
	flag_final = "no"
	for i in range(8,1,-1):
		counts = list(kcounter.count_kmers(seq, i, canonical_kmers=False).values())
		max_c = max(counts)
		ratio = max_c / sum(counts)
		flag = "no"
		if i == 2 and ratio >= 0.3:
			flag="yes"
		elif i == 3 and ratio >= 0.2:
			flag="yes"
		elif i == 4 and ratio >= 0.15:
			flag="yes"
		elif i == 5 and max_c >= 3:
			flag="yes"
		elif i == 6 and max_c >= 2:
			flag="yes"
		elif i == 7 and max_c >= 2:
			flag="yes"
		elif i == 8 and max_c >= 2:
			flag="yes"
		print(f"{seq} // {i} // {max(counts)} // {sum(counts)} // {ratio} // {counts}")
		out.writerow([seq,i,max(counts),sum(counts),ratio,flag,counts])
		# print(f"{seq}\t{i}\t{max(counts)}\t{sum(counts)}\t{ratio}\t{counts}")
		# if ratio > 0.3:
		if flag == "yes" and flag_final == "no":
			flag_final = "low_complexity_kmer" + str(i)
			# break
	return flag_final


def set_default(args):
	args["input_repeats"] = "/clusterfs/jgi/groups/science/metagen/virus/database/spacerextractor_db/Db_SE_GTDBr214_GB_IMG/Repeats.fna"
	args["ref_file"] = "/clusterfs/jgi/groups/science/metagen/virus/database/spacerextractor_db/Db_SE_GTDBr214_GB_IMG/Arrays_info.tsv"
	args["out_file"] = "Repeats_complexity.tsv"
	args["out_file_detail"] = "Repeats_kmer_complexity_detailed.tsv"


if __name__ == "__main__":
	output = main()
