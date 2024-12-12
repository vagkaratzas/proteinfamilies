#!/usr/bin/env python

import sys
import argparse
import os
import gzip
from Bio import AlignIO, SeqIO

def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-a",
        "--alignment",
        required=True,
        metavar="FILE",
        type=str,
        help="Stockholm format multiple sequence alignment from hmmsearch.",
    )
    parser.add_argument(
        "-d",
        "--domtbl",
        required=True,
        metavar="FILE",
        type=str,
        help="Domain summary annotations result from hmmsearch.",
    )
    parser.add_argument(
        "-l",
        "--length_threshold",
        required=True,
        metavar="FLOAT",
        type=float,
        help="Minimum length percentage threshold of annotated domain (env) against query to keep.",
    )
    parser.add_argument(
        "-m",
        "--out_msa",
        required=True,
        metavar="FILE",
        type=str,
        help="Name of the output fasta file with the fasta converted multiple sequence alignment (includes gaps).",
    )
    parser.add_argument(
        "-o",
        "--out_fasta",
        required=True,
        metavar="FILE",
        type=str,
        help="Name of the output fasta file with the fasta converted sequences (no gaps).",
    )
    return parser.parse_args(args)

def filter_sequences(domtbl, length_threshold):
    filtered_sequences = []
    with gzip.open(domtbl, 'rt', encoding='utf-8') as file:
        for line in file:
            if line.startswith('#'):
                continue  # Skip comments

            columns = line.split()
            try:
                qlen = float(columns[5])
                env_from = int(columns[19])
                env_to = int(columns[20])
                env_length = env_to - env_from + 1

                if env_length >= length_threshold * qlen:
                    sequence_name = columns[0]
                    ali_from = int(columns[17])
                    ali_to = int(columns[18])

                    filtered_sequences.append(f"{sequence_name}/{ali_from}-{ali_to}")
            except (IndexError, ValueError):
                continue  # Skip malformed lines

    return filtered_sequences

def filter_stockholm_to_fasta(alignment, filtered_sequences, out_msa, out_fasta):
    base_filename = os.path.basename(alignment).split('.')[0]

    with gzip.open(alignment, 'rt', encoding='utf-8') as file:
        alignment_data = AlignIO.read(file, "stockholm")
        filtered_records = [record for record in alignment_data if record.id in filtered_sequences]

        for record in filtered_records:
            ungapped_length = len(record.seq.replace('-', ''))
            # Extract the part after the slash in the name (if present)
            if '/' in record.id:
                name, range_info = record.id.split('/')
                start, end = map(int, range_info.split('-'))
                # If the sequence starts at 1 and the end matches the ungapped length, remove the slash and range
                if start == 1 and end == ungapped_length:
                    record.id = name

            # Add the family name to the description field
            record.description = base_filename

        with gzip.open(out_msa, 'wt') as gz_file:
            SeqIO.write(filtered_records, gz_file, "fasta")

        for record in filtered_records:
            record.letter_annotations = {}  # Clear annotations
            record.seq = record.seq.replace('-', '')

        with gzip.open(out_fasta, 'wt') as gz_file:
            SeqIO.write(filtered_records, gz_file, "fasta")

        print(f"Filtered alignment saved to {out_msa}")

def filter_recruited(alignment, domtbl, length_threshold, out_msa, out_fasta):
    filtered_sequences = filter_sequences(domtbl, length_threshold)
    filter_stockholm_to_fasta(alignment, filtered_sequences, out_msa, out_fasta)

def main(args=None):
    args = parse_args(args)
    filter_recruited(args.alignment, args.domtbl, args.length_threshold, args.out_msa, args.out_fasta)

if __name__ == "__main__":
    sys.exit(main())
