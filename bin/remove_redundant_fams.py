#!/usr/bin/env python

import sys
import gzip
import argparse
import csv
import os
import shutil
from Bio import SeqIO

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

def remove_redundant_fams(mapping, domtbl, fasta_folder, length_threshold, out_folder):
    for file_name in os.listdir(fasta_folder):
        source_file = os.path.join(fasta_folder, file_name)
        destination_file = os.path.join(out_folder, file_name)

        # TODO logic with actual family filtering
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
