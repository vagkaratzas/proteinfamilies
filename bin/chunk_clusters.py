#!/usr/bin/env python

import sys
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
    # TODO


if __name__ == "__main__":
    sys.exit(main())
