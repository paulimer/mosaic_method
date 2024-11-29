#!/usr/bin/env python3
import argparse
import concurrent.futures
import os
import re
import subprocess as sp
import tempfile

import Bio.Seq
import Bio.SeqIO
import pandas as pd


def remove_AmbiguousIUPAC(genome_path, output_dir):
    """Replaces letters RYWSMKHBVD by N"""
    records = []
    with open(genome_path, 'r') as genome_handle:
        for record in Bio.SeqIO.parse(genome_handle, 'fasta'):
            seq = re.sub(r'[RYWSMKHBVD]', 'N', str(record.seq))
            record.seq = Bio.Seq.Seq(seq)
            records.append(record)
    output_path = os.path.join(output_dir, os.path.basename(genome_path))
    with open(output_path, 'w') as output_handle:
        Bio.SeqIO.write(records, output_handle, 'fasta')


def mask_repeats(genome_path, output_dir, lastz_tools_dir):
    """Masks repeats in a genome"""
    genome_name = os.path.basename(genome_path)
    output_dir = output_dir.rstrip("/")
    output_path = os.path.join(output_dir, genome_name)
    masked_intervals_path = os.path.join(output_dir, genome_name + ".masked_intervals.dat")
    fasta_fragments_py = os.path.join(lastz_tools_dir, "fasta_fragments.py")
    fasta_softmask_intervals_py = os.path.join(lastz_tools_dir, "fasta_softmask_intervals.py")
    detect_repeats_cmd = f"cat {genome_path} | python {fasta_fragments_py} --fragment=200 --step=100 | lastz {genome_path}[multiple,unmask,nameparse=darkspace] /dev/stdin --masking=3 --progress+masking=10K --format=none --outputmasking+:soft={masked_intervals_path}"
    mask_repeats_cmd = f"cat  {genome_path} | python {fasta_softmask_intervals_py} --origin=1 {masked_intervals_path} > {output_path}"
    sp.run(detect_repeats_cmd, shell=True, check=True)
    sp.run(mask_repeats_cmd, shell=True, check=True)



def pretreat_genomes(directory, taxon_csv, output, lastz_tools_dir, threads=1):
    """Masks repeats and replaces ambiguous IUPAC letters by N"""
    taxon_df = pd.read_csv(taxon_csv)
    genomes = taxon_df["genome"].tolist()
    genomes = [os.path.join(directory, g) for g in genomes if os.path.exists(os.path.join(directory, g))]
    n_genomes = len(genomes)
    os.makedirs(output, exist_ok=True)
    with tempfile.TemporaryDirectory() as temp_dir:
        original_genomes = [os.path.join(directory, g) for g in genomes]
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
            executor.map(remove_AmbiguousIUPAC, original_genomes, [temp_dir]*n_genomes)
        temp_genomes = [os.path.join(temp_dir, g) for g in genomes]
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
            executor.map(mask_repeats, temp_genomes, [output]*n_genomes, [lastz_tools_dir]*n_genomes)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Masks repeats and replaces ambiguous IUPAC letters by N")
    parser.add_argument(
        "--threads",
        "-t",
        type=int,
        default=1,
        help="Number of threads to use"
    )
    parser.add_argument(
        "lastz_tools_dir",
        type=str,
        help="Directory containing lastz tools"
    )
    parser.add_argument(
        "directory",
        type=str,
        help="Directory containing genomes to pretreat"
    )
    parser.add_argument(
        "taxon_csv",
        type=str,
        help="CSV file containing genomes and their taxon"
    )
    parser.add_argument(
        "output",
        type=str,
        help="Output directory"
    )
    args = parser.parse_args()
    pretreat_genomes(args.directory, args.taxon_csv, args.output, args.lastz_tools_dir, args.threads)
