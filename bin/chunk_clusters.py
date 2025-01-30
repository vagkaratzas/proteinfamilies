#!/usr/bin/env python

import sys
import os
import argparse
import polars as pl
from Bio import SeqIO


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--clustering", required=True, type=str, help="TSV clustering file input.")
    parser.add_argument("-s", "--sequences" , required=True, type=str, help="Initial sequences FASTA file.")
    parser.add_argument("-t", "--threshold" , required=True, type=int, help="Minimum cluster size to keep.")
    parser.add_argument("-o", "--out_folder", required=True, type=str, help="Output folder to store FASTA files.")
    return parser.parse_args(args)


def collect_clusters(clustering_file, threshold):
    """Load clusters efficiently using Polars and filter by threshold."""
    df = pl.read_csv(
        clustering_file,
        separator="\t",
        has_header=False,
        new_columns=["rep", "member"],
        schema_overrides={"rep": pl.Utf8, "member": pl.Utf8},  # Treat both columns as strings
    )

    # Group by representative, collect members, filter by threshold
    clusters = (
        df.group_by("rep")
        .agg(pl.col("member"))
        .filter(pl.col("member").list.len() >= threshold)
    )

    # Convert to dictionary {rep: frozenset(members)}
    return {row["rep"]: frozenset(row["member"]) for row in clusters.iter_rows(named=True)}


def process_clusters(clusters, fasta_file, out_folder):
    """Write FASTA sequences efficiently using an indexed approach."""
    os.makedirs(out_folder, exist_ok=True)

    # Use SeqIO.index() to create an indexed FASTA dictionary
    seq_dict = SeqIO.index(fasta_file, "fasta")

    for chunk_num, (rep, members) in enumerate(clusters.items(), start=1):
        output_file = os.path.join(out_folder, f"{chunk_num}.fasta")

        # Collect sequences that belong to this cluster
        cluster_seqs = [seq_dict[member] for member in members if member in seq_dict]

        if cluster_seqs:
            with open(output_file, "w") as fasta_out:
                SeqIO.write(cluster_seqs, fasta_out, "fasta")


def main(args=None):
    args = parse_args(args)

    # Collect valid clusters with Polars
    clusters = collect_clusters(args.clustering, args.threshold)

    # Process and write FASTA sequences
    process_clusters(clusters, args.sequences, args.out_folder)


if __name__ == "__main__":
    sys.exit(main())
