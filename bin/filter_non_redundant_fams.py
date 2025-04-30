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
        "-r",
        "--redundant_ids",
        required=True,
        metavar="FILE",
        type=str,
        help="Text file with one redundant family ID per line.",
    )
    parser.add_argument(
        "-s",
        "--seqs",
        required=True,
        metavar="FOLDER",
        type=str,
        help="All family fasta sequence files.",
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
        "-se",
        "--seeds",
        required=True,
        metavar="FOLDER",
        type=str,
        help="All family seed MSAs.",
    )
    parser.add_argument(
        "-a",
        "--alns",
        required=True,
        metavar="FOLDER",
        type=str,
        help="All family full MSAs.",
    )
    return parser.parse_args(args)


def read_redundant_ids(filepath):
    with open(filepath) as f:
        return set(line.strip() for line in f if line.strip())


def filter_files(input_dir, output_dir, redundant_ids):
    os.makedirs(output_dir, exist_ok=True)
    for file in os.listdir(input_dir):
        fam_id = file.split(".")[0]
        if fam_id not in redundant_ids:
            shutil.copy(
                os.path.join(input_dir, file),
                os.path.join(output_dir, file)
            )


def main(args=None):
    args = parse_args(args)
    redundant_ids = read_redundant_ids(args.redundant_ids)

    filter_files(args.seqs, "fasta", redundant_ids)
    filter_files(args.models, "hmm", redundant_ids)
    filter_files(args.seeds, "seed_msa", redundant_ids)
    filter_files(args.alns, "full_msa", redundant_ids)


if __name__ == "__main__":
    sys.exit(main())
