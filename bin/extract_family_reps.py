#!/usr/bin/env python

import sys
import os
import gzip
import argparse
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
        "-o",
        "--out_file",
        required=True,
        metavar="FILE",
        type=str,
        help="Name of the output fasta file with family representative sequences.",
    )
    return parser.parse_args(args)

def extract_first_sequences(msa_folder, out_file):
    # Open the output file in write mode
    with open(out_file, "w") as outfile:
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
                    # Create a new SeqRecord with the cleaned sequence
                    cleaned_record = SeqRecord(
                        Seq(cleaned_sequence),
                        id=first_record.id,
                        description=first_record.description
                    )
                    # Write the cleaned sequence to the output file
                    SeqIO.write(cleaned_record, outfile, "fasta")

def main(args=None):
    args = parse_args(args)
    extract_first_sequences(args.full_msa_folder, args.out_file)

if __name__ == "__main__":
    sys.exit(main())
