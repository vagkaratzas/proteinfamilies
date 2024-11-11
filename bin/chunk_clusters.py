#!/usr/bin/env python

import sys
import os
import argparse
import pandas as pd


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


def main(args=None):
    args = parse_args(args)

    # Load the clustering file
    clustering_df = pd.read_csv(args.clustering, sep='\t', header=None, names=['representative', 'member'])

    # Group by cluster representative and filter by threshold
    grouped = clustering_df.groupby('representative')
    clusters = {rep: members for rep, members in grouped if len(members) >= args.threshold}

    # Create the output folder if it doesn't exist
    os.makedirs(args.out_folder, exist_ok=True)

    # Save each cluster to a separate TSV file
    for i, (rep, members) in enumerate(clusters.items(), start=1):
        output_file = os.path.join(args.out_folder, f"{i}.tsv")
        members.to_csv(output_file, sep='\t', index=False, header=False)
        print(f"Saved cluster {rep} with {len(members)} members to {output_file}")


if __name__ == "__main__":
    sys.exit(main())
