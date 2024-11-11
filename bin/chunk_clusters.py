#!/usr/bin/env python

import sys
import os
import argparse
import pandas as pd
from Bio import SeqIO

def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c",
        "--clustering",
        required=True,
        metavar="FILE",
        type=str,
        help="TSV clustering file input.",
    )
    parser.add_argument(
        "-s",
        "--sequences",
        required=True,
        metavar="FILE",
        type=str,
        help="Initial sequences fasta file.",
    )
    parser.add_argument(
        "-t",
        "--threshold",
        required=True,
        metavar="INT",
        type=int,
        help="Minimum cluster size to keep.",
    )
    parser.add_argument(
        "-o",
        "--out_folder",
        required=True,
        metavar="FOLDER",
        type=str,
        help="Name of the output folder to be created.",
    )
    return parser.parse_args(args)

def load_sequences(sequences_file):
    # Load sequences from the input FASTA file into a dictionary for quick lookup
    return {record.id: record for record in SeqIO.parse(sequences_file, "fasta")}

def main(args=None):
    args = parse_args(args)

    # Create output directory if it doesn't exist
    os.makedirs(args.out_folder, exist_ok=True)

    # Load the sequences from the input FASTA file
    sequences = load_sequences(args.sequences)

    # Read clustering file and filter by threshold
    clusters = pd.read_csv(args.clustering, sep='\t', header=None, names=["rep", "member"])
    cluster_groups = clusters.groupby("rep")

    # Process each cluster and output sequences to chunked FASTA files
    chunk_num = 1
    for rep, group in cluster_groups:
        members = group["member"].tolist()

        # Filter by threshold
        if len(members) >= args.threshold:
            output_file = os.path.join(args.out_folder, f"{chunk_num}.fasta")

            # Write sequences to the chunked FASTA file
            with open(output_file, "w") as fasta_out:
                for member in members:
                    if member in sequences:
                        SeqIO.write(sequences[member], fasta_out, "fasta")
                    else:
                        print(f"Warning: Sequence {member} not found in the input FASTA file.")

            chunk_num += 1

if __name__ == "__main__":
    sys.exit(main())
