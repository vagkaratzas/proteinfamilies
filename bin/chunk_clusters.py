#!/usr/bin/env python

import sys
import os
import argparse
import csv
from collections import defaultdict
from Bio import SeqIO


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--clustering", required=True, type=str, help="TSV clustering file input.")
    parser.add_argument("-s", "--sequences" , required=True, type=str, help="Initial sequences FASTA file.")
    parser.add_argument("-t", "--threshold" , required=True, type=int, help="Minimum cluster size to keep.")
    parser.add_argument("-o", "--out_folder", required=True, type=str, help="Output folder to store FASTA files.")
    return parser.parse_args(args)


def collect_clusters(clustering_file, threshold):
    """Load clusters from file and filter by threshold."""
    clusters = defaultdict(set)  # Use sets for fast lookups

    with open(clustering_file) as f:
        csv_reader = csv.reader(f, delimiter="\t")
        for rep, member in csv_reader:
            clusters[rep].add(member)

    # Filter clusters by threshold
    return {rep: members for rep, members in clusters.items() if len(members) >= threshold}


def load_sequences(fasta_file):
    """Load all sequences into a dictionary (ID -> SeqRecord)."""
    return {record.id: record for record in SeqIO.parse(fasta_file, "fasta")}


def write_fasta(output_file, records):
    """Write multiple FASTA records efficiently."""
    with open(output_file, "w") as fasta_out:
        SeqIO.write(records, fasta_out, "fasta")


def main(args=None):
    args = parse_args(args)

    # Ensure output directory exists
    os.makedirs(args.out_folder, exist_ok=True)

    # Collect valid clusters
    clusters = collect_clusters(args.clustering, args.threshold)

    # Load all sequences into memory for fast lookups
    seq_dict = load_sequences(args.sequences)

    # Process clusters in one pass
    for chunk_num, (rep, members) in enumerate(clusters.items(), start=1):
        output_file = os.path.join(args.out_folder, f"{chunk_num}.fasta")

        # Collect sequences that belong to this cluster
        cluster_seqs = [seq_dict[member] for member in members if member in seq_dict]

        if cluster_seqs:
            write_fasta(output_file, cluster_seqs)


if __name__ == "__main__":
    sys.exit(main())
