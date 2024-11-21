#!/usr/bin/env python

import sys
import os
import gzip
import argparse
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f",
        "--full_msa_folder",
        required=True,
        metavar="FOLDER",
        type=str,
        help="Input folder with Stockholm full alignments.",
    )
    parser.add_argument(
        "-m",
        "--map",
        required=True,
        metavar="FILE",
        type=str,
        help="Name of the output csv file with representative sequences to family id mapping.",
    )
    parser.add_argument(
        "-o",
        "--out_fasta",
        required=True,
        metavar="FILE",
        type=str,
        help="Name of the output fasta file with family representative sequences.",
    )
    return parser.parse_args(args)

def extract_first_sequences(msa_folder, map_file, out_fasta):
    # Open the output FASTA file in write mode
    with open(out_fasta, "w") as fasta_out, open(map_file, "w", newline="") as csv_out:
        csv_writer = csv.writer(csv_out)
        # Write the CSV header
        csv_writer.writerow(["family", "representative_sequence"])

        # Iterate over all files in the MSA folder
        for filename in os.listdir(msa_folder):
            filepath = os.path.join(msa_folder, filename)
            # Parse the Stockholm file and extract the first sequence
            with gzip.open(filepath, "rt") as sto_file:
                records = list(SeqIO.parse(sto_file, "stockholm"))
                if records:
                    first_record = records[0]
                    # Remove gaps from the sequence
                    cleaned_sequence = str(first_record.seq).replace("-", "").replace(".", "")
                    # Modify the ID to only include the part before the first space
                    cleaned_id = first_record.id.split(" ")[0]
                    # Create a new SeqRecord with the cleaned sequence and ID
                    cleaned_record = SeqRecord(
                        Seq(cleaned_sequence),
                        id=cleaned_id,
                        description=""
                    )
                    # Write the cleaned sequence to the FASTA file
                    SeqIO.write(cleaned_record, fasta_out, "fasta")
                    # Write the mapping to the CSV file
                    csv_writer.writerow([os.path.splitext(os.path.splitext(filename)[0])[0], cleaned_id])

def main(args=None):
    args = parse_args(args)
    extract_first_sequences(args.full_msa_folder, args.map, args.out_fasta)

if __name__ == "__main__":
    sys.exit(main())
