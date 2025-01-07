#!/usr/bin/env python

import sys
import numpy as np
import argparse
from Bio import SeqIO
from Bio.Seq import Seq


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-a",
        "--alignment",
        required=True,
        metavar="FILE",
        type=str,
        help="Multiple sequence alignment file in fasta format.",
    )
    parser.add_argument(
        "-g",
        "--gap_threshold",
        required=True,
        metavar="FLOAT",
        type=float,
        help="Minimum gap occupancy across sequences to keep.",
    )
    parser.add_argument(
        "-o",
        "--out_fasta",
        required=True,
        metavar="FILE",
        type=str,
        help="Name of the output fasta file with the trimmed multiple sequence alignment.",
    )
    return parser.parse_args(args)


def read_fasta_to_matrix(file_path):
    records = list(SeqIO.parse(file_path, "fasta"))
    max_length = max(len(record.seq) for record in records)
    matrix = np.zeros((len(records), max_length), dtype=np.dtype("U1"))
    original_names = []

    for i, record in enumerate(records):
        original_names.append(record.id)
        matrix[i, : len(record.seq)] = list(str(record.seq))

    return matrix, original_names


def calculate_trim_positions(sequence_matrix, gap_threshold):
    numeric_matrix = np.where(sequence_matrix == "-", 0, 1)
    num_rows = numeric_matrix.shape[0]
    column_sums = np.sum(numeric_matrix, axis=0)
    column_sums_percentage = column_sums / num_rows
    start_position = np.argmax(column_sums_percentage > gap_threshold)
    end_position = (
        len(column_sums_percentage)
        - np.argmax(column_sums_percentage[::-1] > gap_threshold)
        - 1
    )

    return start_position, end_position


def write_trimmed_sequences(
    sequence_matrix_trimmed, original_sequence_names, out_fasta
):
    trimmed_records = []
    for i, sequence in enumerate(sequence_matrix_trimmed):
        trimmed_sequence = "".join(map(str, sequence))
        original_name = original_sequence_names[i]
        trimmed_record = SeqIO.SeqRecord(
            Seq(trimmed_sequence), id=original_name, description=""
        )
        trimmed_records.append(trimmed_record)

    with open(out_fasta, "w") as output_fasta:
        SeqIO.write(trimmed_records, output_fasta, "fasta")


def trim_msa(alignment, gap_threshold, out_fasta):
    sequence_matrix, original_sequence_names = read_fasta_to_matrix(alignment)
    start_position, end_position = calculate_trim_positions(
        sequence_matrix, gap_threshold
    )
    sequence_matrix_trimmed = sequence_matrix[:, start_position : end_position + 1]
    write_trimmed_sequences(sequence_matrix_trimmed, original_sequence_names, out_fasta)


def main(args=None):
    args = parse_args(args)
    trim_msa(args.alignment, args.gap_threshold, args.out_fasta)


if __name__ == "__main__":
    sys.exit(main())
