#!/usr/bin/env python
import sys
import argparse
import gzip
from Bio import SeqIO, SearchIO


def open_file(filename):
    """Helper function to open both gzipped and regular files"""
    if filename.endswith(".gz"):
        return gzip.open(filename, "rt")
    return open(filename, "r")


def main():
    parser = argparse.ArgumentParser(
        description="Get the sequences that are missing from the hmmscan txt output table."
    )
    parser.add_argument(
        "-d",
        "--hmmscantxt",
        required=True,
        metavar="FILE",
        type=str,
        help="The hmmsearch txt output (gzipped).",
    )
    parser.add_argument(
        "-f",
        "--fasta",
        required=True,
        metavar="FILE",
        type=str,
        help="The FASTA (gzip or not) file containing all the sequences.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        metavar="FILE",
        type=str,
        help="Output gzipped FASTA file path for sequences without matches.",
    )
    args = parser.parse_args()

    # Collect sequence IDs with matches
    hit_sequence_ids = set()
    with gzip.open(args.hmmscantxt, "rt") as domtable_fh:
        for hmmer_result in SearchIO.parse(domtable_fh, "hmmer3-text"):
            hit_sequence_ids.add(hmmer_result.id)

    with open_file(args.fasta) as fasta_fl, gzip.open(
        args.output, "wt"
    ) as fasta_out_fl:

        output_sequences = []
        for record in SeqIO.parse(fasta_fl, "fasta"):
            if record.id not in hit_sequence_ids:
                output_sequences = [record]

        SeqIO.write(output_sequences, fasta_out_fl, "fasta")


if __name__ == "__main__":
    sys.exit(main())
