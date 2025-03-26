#!/usr/bin/env python

## Originally written by Evangelos Karatzas and released under the MIT license.
## See git repository (https://github.com/nf-core/proteinfamilies) for full license text.

import os
import sys
import argparse
import shutil


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-s",
        "--seqs",
        required=True,
        metavar="FOLDER",
        type=str,
        help="Filtered fasta files to grab names from.",
    )
    parser.add_argument(
        "-m",
        "--models",
        required=True,
        metavar="FOLDER",
        type=str,
        help="All family HMMs.",
    )
    parser.add_argument(
        "-o",
        "--out_folder",
        required=True,
        metavar="FOLDER",
        type=str,
        help="Name of the output folder file with the filtered HMMs.",
    )
    return parser.parse_args(args)


def filter_non_redundant_hmms(seqs, models, out_folder):
    seq_basenames = {
        os.path.basename(f).split(".")[0]
        for f in os.listdir(seqs)
        if f.endswith(".fasta.gz") or f.endswith(".fasta")
    }
    # Iterate through the models folder and copy matching files
    for model_file in os.listdir(models):
        if model_file.endswith(".hmm.gz"):
            model_basename = os.path.basename(model_file).split(".")[0]
            model_chunk = model_basename.split("_")[-1] # for the case where we didn't fish additional sequences, fastas only have the chunk id for name
            if model_basename in seq_basenames or model_chunk in seq_basenames:
                src = os.path.join(models, model_file)
                dst = os.path.join(out_folder, model_file)
                shutil.copy(src, dst)


def main(args=None):
    args = parse_args(args)

    os.makedirs(args.out_folder, exist_ok=True)
    filter_non_redundant_hmms(args.seqs, args.models, args.out_folder)


if __name__ == "__main__":
    sys.exit(main())
