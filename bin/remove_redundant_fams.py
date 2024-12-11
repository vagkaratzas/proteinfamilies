#!/usr/bin/env python

import os
import sys
import pandas as pd
import argparse
import shutil

def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-m",
        "--mapping",
        required=True,
        metavar="FILE",
        type=str,
        help="CSV metadata mapping input.",
    )
    parser.add_argument(
        "-d",
        "--domtbl",
        required=True,
        metavar="FILE",
        type=str,
        help="TSV hmmsearch domtbl out results for filtering.",
    )
    parser.add_argument(
        "-f",
        "--fasta_folder",
        required=True,
        metavar="FOLDER",
        type=str,
        help="Name of the input folder file with the pre-filtered fasta.",
    )
    parser.add_argument(
        "-l",
        "--length_threshold",
        required=True,
        metavar="FLOAT",
        type=str,
        help="Minimum length percentage threshold of annotated domain (env) against query to keep.",
    )
    parser.add_argument(
        "-o",
        "--out_folder",
        required=True,
        metavar="FOLDER",
        type=str,
        help="Name of the output folder file with the filtered fasta.",
    )
    return parser.parse_args(args)

def remove_self_hits(mapping_set, domtbl_df):
    # Filter out these self-hits from domtbl_df based on the set membership
    filtered_domtbl_df = domtbl_df[
        ~domtbl_df.apply(lambda row: (row["target name"], row["query name"]) in mapping_set, axis=1)
    ]

    return filtered_domtbl_df

def remove_redundant_fams(mapping, domtbl, fasta_folder, length_threshold, out_folder):
    mapping_df = pd.read_csv(
        mapping,
        comment='#',
        usecols=["Family Id", "Size", "Representative Id"]
    )
    domtbl_df = pd.read_csv(
        domtbl,
        sep=r'\s+',
        comment='#',
        header=None,
        usecols=[0, 3, 5, 20, 21]
    ).rename(columns={0: "target name", 3: "query name", 5: "qlen", 20: "env from", 21: "env to"})
    # Create a set of (Representative Id, Family Id) pairs from mapping_df
    mapping_set = set(zip(mapping_df["Representative Id"], mapping_df["Family Id"]))

    domtbl_df = remove_self_hits(mapping_set, domtbl_df)
    
    print(domtbl_df)

    # TODO logic with actual family filtering
    # 1. remove self-hits, done
    # 2. filter_by_length
    # 3. keep larger hit
    # 4. write out non redundant
    for file_name in os.listdir(fasta_folder):
        source_file = os.path.join(fasta_folder, file_name)
        destination_file = os.path.join(out_folder, file_name)

        # Check if it is a file (not a directory)
        if os.path.isfile(source_file):
            # Copy file to the destination folder
            shutil.copy2(source_file, destination_file)

def main(args=None):
    args = parse_args(args)

    os.makedirs(args.out_folder, exist_ok=True)
    remove_redundant_fams(args.mapping, args.domtbl, args.fasta_folder, args.length_threshold, args.out_folder)

if __name__ == "__main__":
    sys.exit(main())
