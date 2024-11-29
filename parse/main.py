#!/usr/bin/env python3

import argparse
import sys

import numpy as np
import pandas as pd
from fun import *

def main():
    parser = argparse.ArgumentParser(
        description="""
        Concatenates single mlds from alignements, saves them. Optionnaly.
        """)
    parser.add_argument(
        "--pairs",
        action="store_true",
        help="only align pairs of genes from coresimul simulation"
    )
    parser.add_argument(
        "--threads",
        type=int,
        help="Number of threads to use for parsing csv files",
        default=1
    )
    parser.add_argument(
        "--from_sqlite_db",
        type=str,
        help="Start from sqlite database. Incompatible with from_florian_mld or from_csv"
    )
    parser.add_argument(
        "--taxon_csv",
        help="File defining genome-taxon relationships",
        type=str
    )
    parser.add_argument(
        "--cluster_name",
        help="The column name (genus, family, or custom cluster name) of the taxon csv file to consider for comparisons",
        type=str
    )
    parser.add_argument(
        "--full_mld",
        help="the full mld directory",
        type=str
    )
    parser.add_argument(
        "--binned_mld",
        help="the binned mld directory",
        type=str
    )
    args = parser.parse_args()

    if args.from_sqlite_db:
        taxon_df = pd.read_csv(args.taxon_csv)
        level_list = sorted(taxon_df[args.cluster_name].unique())
        levels = list(itertools.combinations(level_list, 2))
        for level in levels:
            genome_comps = get_genome_comp(level, args.taxon_csv, "", args.cluster_name, output_csv=False, pairs=args.pairs)
            full_mld = get_all_mlds(genome_comps, args.from_sqlite_db, threads=args.threads)
            summed_mld = sum_mlds(full_mld)
            binned_mld = bin_mld(
                summed_df=summed_mld,
                linear_bin_width=3,
                limit_size=35.5,
                power_increment=0.1,
                ncomp=full_mld.shape[0]
            )

            full_mld.to_csv(os.path.join(args.full_mld, f"{level[0]}_vs_{level[1]}_full_mld_comp.csv"), index=False)
            binned_mld.to_csv(os.path.join(args.binned_mld, f"{level[0]}_vs_{level[1]}_binned_mld.csv"), index=False)


if __name__ == "__main__":
    main()
