#!/usr/bin/env python
import sys
import csv
import os
import argparse
from Bio.Seq import Seq


def main():
    print("starting")
    n = 0
    cluster = {}
    print("reading tsv")
    parser = argparse.ArgumentParser()
    parser.set_defaults(func=main)
    parser.set_defaults(program="predict")
    parser.add_argument('--in_file','-i', dest='in_file_spacers', required=False, default='/tmp/All_spacers_info_filtered-Jul19-24.tsv',help='input tsv file of spacers')
    parser.add_argument('--out_file','-o', dest='clstr_file', required=False, default='/tmp/All_spacers_info_filtered_clusters-Jul19-24.tsv',help='output file for clusters')
    args = vars(parser.parse_args())
    
    in_file_spacers = args["in_file_spacers"]
    clstr_file = args["clstr_file"]

    with open(in_file_spacers, "r") as tsv_file, open(clstr_file, "w") as out_file:
        tsv_file = csv.reader(tsv_file, delimiter="\t")
        tsv_writer = csv.writer(out_file, delimiter='\t', quotechar=' ', quoting=csv.QUOTE_NONE)
        for line in tsv_file:
            if line[0] == "spacer_id":
                tsv_writer.writerow(["cluster_id", "spacer_id"])
                # print("Spacer id")
                continue
            seq = line[1]
            revseq = Seq(seq).reverse_complement()
            if seq in cluster:
                tsv_writer.writerow([cluster[seq], line[0]])
                # cluster[seq].append(line[0])
            elif revseq in cluster:
                tsv_writer.writerow([cluster[revseq], line[0]])
                # cluster[revseq].append(line[0])
            else:
                n = n + 1
                cluster_id = f"Sp_cl_{n:012d}"
                if n % 100000 == 0:
                    print(f"new cluster .. {cluster_id}")
                cluster[seq] = cluster_id
                # print(f"{cluster_id} \\ {line[0]}")
                tsv_writer.writerow([cluster_id, line[0]])
    print("finished")



if __name__ == "__main__":
    output = main()
