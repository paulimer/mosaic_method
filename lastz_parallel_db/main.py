#!/usr/bin/env python3

import argparse

from utils import *

def main():
    parser = argparse.ArgumentParser(
        description="""
        Runs lastz, saves its output as smaller csv files.
        """
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=10,
        help="Number of threads to use (default: 10)"
    )
    parser.add_argument(
        "-u",
        "--update",
        action="store_true",
        help="Update the existing database instead of creating a new one"
    )
    parser.add_argument(
        "-p",
        "--pairs",
        action="store_true",
        help="only align pairs of genes from coresimul simulation"
    )
    parser.add_argument(
        "-a",
        "--aligner",
        type=str,
        default="lastz",
        help="Aligner to use (default: lastz)"
    )
    parser.add_argument(
        "taxon_csv",
        type=str,
        help="Path to the taxon csv file"
    )
    parser.add_argument(
        "genomes_path",
        type=str,
        help="Path to the genomes directory"
    )
    parser.add_argument(
        "cluster_name",
        type=str,
        help="The name of the cluster level(column in taxon_csv) to use to make comparisons"
    )
    parser.add_argument(
        "output_db",
        type=str,
        help="The output sqlite file"
    )
    args = parser.parse_args()
    if not args.aligner in ["lastz", "mummer"]:
        print("Invalid aligner")
        return
    con = create_lastz_db(args.taxon_csv, args.genomes_path, args.cluster_name, args.output_db, args.threads, args.update, args.pairs, args.aligner)
    con.close()

if __name__ == "__main__":
    main()
