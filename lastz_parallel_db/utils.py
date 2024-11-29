import concurrent.futures
import io
import itertools
import os
import re
import subprocess as sp
import sys
import tempfile

import sqlite3
import pandas as pd
import numpy as np


def parse_cigarx_line(line):
    """
    Parses a cigarx line and counts the length of matches.
    """
    li = len(line) - 1
    matches = []
    while li >= 0:
        if line[li] == "=":
            len_match = ""
            li -= 1
            while line[li].isdigit():
                len_match = line[li] + len_match
                li -= 1
                if li == -1:
                    break
            if len_match:
                matches.append(int(len_match))
            else:
                matches.append(1)
        else:
            li -= 1
    count_array = np.zeros(max(matches), dtype=int)
    for match_length in matches:
        # first element is matches of length 1
        count_array[match_length - 1] += 1
    return count_array


def conc_show_aligns_matches(mummer_out):
    """
    Concatenates all matches and mismatches lines from a mummer show-aligns alignment.
    """

    # concatenate all lines
    misandmatches = ""

    # is the previous line an alignment line (with atgc)
    prev_align = False
    cur_align = False

    # calc the padding before the alignment
    counted_start_letter = False
    padding = 0


    lines = mummer_out.splitlines()
    for line in lines:
        if line == "":
            continue
        prev_align = cur_align
        if line[0].isdigit():
            cur_align = True
        else:
            cur_align = False
        if not counted_start_letter and cur_align:
            counted_start_letter = True
            detect_bases = re.compile(r"[atgcATGC]")
            base_match = detect_bases.search(line)
            if base_match:
                padding = base_match.start()

        if prev_align and not cur_align:
            pattern = re.compile(r"[^ \^]")
            if pattern.search(line):
                print("Error in line : ")
                print(line)
                break
            misandmatches += line[padding:]
    return misandmatches


def mummer_count_matches(misandmatches):
    """
    Counts the length of matches in a string of matches and mismatches.
    """
    matches = []
    cur_match = 0
    for char in misandmatches:
        if char == "^":
            if cur_match:
                matches.append(cur_match)
                cur_match = 0
        elif char == " ":
            cur_match += 1
        else:
            print(f"Error with char {char}")
    if cur_match:
        matches.append(cur_match)
    count_array = np.zeros(max(matches), dtype=int)
    for match_length in matches:
        # first element is matches of length 1
        count_array[match_length - 1] += 1
    return count_array


def calc_count_array_len(summed_count_array):
    total_match = 0
    for i, count in enumerate(summed_count_array):
        total_match += (i+1)*count
    return total_match


def parse_show_align_output(mummer_out):
    """
    Splits a mummer show-aligns output in multiple alignments.
    """
    count_arrays = []
    length1s = []
    sum_matches = []
    alignment = ""
    split_align_pattern = r'(?s)(BEGIN.*?\])(.*?)(END.*?)(?=\n|$)'
    divs = []
    for m in re.finditer(split_align_pattern, mummer_out):
        alignment = m.group(2).strip("\n")
        misandmatches = conc_show_aligns_matches(alignment)
        count_array = mummer_count_matches(misandmatches)
        sum_matches.append(calc_count_array_len(count_array))
        count_arrays.append(count_array)
        stats = m.group(1).split(" ")
        length1 = abs(int(stats[6]) - int(stats[4]))
        length1s.append(length1)
        divs.append(calc_count_array_len(count_array)/length1)
    return count_arrays, length1s, sum_matches, max(divs)


def run_mummer(target, query, prefix):
    """
    Runs mummer, adds matches length counts together and returns them as a counter.
    """
    # flag to know if a temporary file was created
    temp_create = False
    sequence_names = re.compile(r">[^\n]*")
    with open(target, "r") as f:
        target_names = sequence_names.findall(f.read())
    with open(query, "r") as f:
        query_names = sequence_names.findall(f.read())
    if len(target_names) > 1:
        temp_create = True
        # concatenate all target sequences in one file
        with open(target, "r") as f:
            target_seq = f.read()
        # remove all sequence names except the first one
        target_seq_lines = target_seq.splitlines()
        target_seq_lines = [target_seq_lines[0]] + [line for line in target_seq_lines[1:] if not line.startswith(">")]
        target_seq = "\n".join(target_seq_lines)
        # tempfile with target sequences concatenated
        target_f = tempfile.NamedTemporaryFile(mode="w", delete=False)
        target_f.write(target_seq)
        target_f.close()
        target = target_f.name

    _ = sp.run(
        ['nucmer', '--mum', "--prefix", prefix, target, query],
        check=True,
        capture_output=True
    )
    # check for empty delta file (no alignment found)
    with open(f"{prefix}.delta", "r") as f:
        if len(f.readlines()) <= 2:
            if temp_create:
                os.remove(target)
            os.remove(f"{prefix}.delta")
            return np.zeros(1, dtype=int), 1, 0, 0, 0
    target_names = [tn[1:] for tn in target_names]
    query_names = [qn[1:] for qn in query_names]
    summed_count_array = None
    summed_matches = 0
    summed_aligned = 0
    max_divs = []
    for query_name in query_names:
        shal_command = ['show-aligns', '-r', f"{prefix}.delta", target_names[0].split(" ")[0], query_name]
        try:
            res_show_aligns = sp.run(
                shal_command,
                capture_output=True,
                check=True,
                encoding="utf-8"
            )
        except sp.CalledProcessError:
            # case where a given contig does not have a single match in the target
            # another possibility would be to scan for them beforehand
            # but better ask for forgiveness
            continue
        count_arrays, length1s, sum_matches, max_div = parse_show_align_output(res_show_aligns.stdout)
        for count_array, length1, sum_match in zip(count_arrays, length1s, sum_matches):
            if summed_count_array is None:
                summed_count_array = count_array
                summed_matches = sum_match
                summed_aligned = length1
                max_divs.append(max_div)
                continue
            if len(count_array) < len(summed_count_array):
                count_array = np.pad(count_array, (0, len(summed_count_array) - len(count_array)), mode='constant')
            elif len(count_array) > len(summed_count_array):
                summed_count_array = np.pad(summed_count_array, (0, len(count_array) - len(summed_count_array)), mode='constant')
            summed_count_array += count_array
            summed_matches += sum_match
            summed_aligned += length1
            max_divs.append(max_div)
    if summed_aligned == 0:
        average_divergence = 1
        summed_count_array = np.zeros(1, dtype=int)
        top_div = 0
    else:
        average_divergence = 1 - summed_matches / summed_aligned
        top_div = max(max_divs)
    if temp_create:
        os.remove(target)
    os.remove(f"{prefix}.delta")
    return summed_count_array, average_divergence, summed_matches, summed_aligned, top_div


def run_lastz(target, query):
    """
    Runs lastz, adds cigarx matches length counts together and returns them as a counter.
    """
    target = target+"[multiple]"
    query = query#+"[multiple]"

    res_lastz = sp.run(
        ['lastz',target,query,"--format=general:length1,idfrac,cigarx"],
        capture_output=True,
        check=True,
        encoding="utf-8"
    )
    res_lastz_df = pd.read_csv(io.StringIO(res_lastz.stdout), sep="\t")
    summed_count_array = None
    summed_matches = 0
    summed_aligned = 0
    divs = []
    for _, row in res_lastz_df.iterrows():
        # cigarx counting
        count_array = parse_cigarx_line(row["cigarx"])
        # Pad the arrays with zeroes if they have different sizes
        if summed_count_array is None:
            summed_count_array = count_array
            summed_matches = int(row["idfrac"].split("/")[0])
            summed_aligned = int(row["idfrac"].split("/")[1])
            divs.append(summed_matches/summed_aligned)
            continue
        if len(count_array) < len(summed_count_array):
            count_array = np.pad(count_array, (0, len(summed_count_array) - len(count_array)), mode='constant')
        elif len(count_array) > len(summed_count_array):
            summed_count_array = np.pad(summed_count_array, (0, len(count_array) - len(summed_count_array)), mode='constant')
        summed_count_array += count_array

        # idfrac calculations
        summed_matches += int(row["idfrac"].split("/")[0])
        summed_aligned += int(row["idfrac"].split("/")[1])
        divs.append(int(row["idfrac"].split("/")[0])/int(row["idfrac"].split("/")[1]))
    if summed_aligned == 0:
        average_divergence = 1
        summed_count_array = np.zeros(1, dtype=int)
        top_div = 0
    else:
        average_divergence = 1 - summed_matches / summed_aligned
        top_div = max(divs)
    return summed_count_array, average_divergence, summed_matches, summed_aligned, top_div


def align_exec(genome_pair, align="lastz", prefix=""):
    """
    computes a lastz entry for two genomes
    """
    genome_1, genome_2 = genome_pair
    if align == "lastz":
        count_array, average_divergence, total_matches, total_aligned, top_div = run_lastz(genome_1, genome_2)
    else:
        count_array, average_divergence, total_matches, total_aligned, top_div = run_mummer(genome_1, genome_2, prefix)
    genome_small = min(genome_1, genome_2)
    genome_large = max(genome_1, genome_2)
    return genome_small, genome_large, count_array, average_divergence, total_matches, total_aligned, top_div


def align_entry(res, con):
    """
    inserts a result in database
    """
    cur = con.cursor()
    res = list(res)
    res[2] = res[2].tobytes()
    res[0] = os.path.basename(res[0])
    res[1] = os.path.basename(res[1])
    cur.execute("INSERT INTO lastz VALUES (?, ?, ?, ?, ?, ?, ?)", res)
    con.commit()


def get_genome_comps(taxon_csv, pairs=False):
    """
    Generates all genome combinations for a list of clusters
    """
    taxon_df = pd.read_csv(taxon_csv)
    if not pairs:
        sorted_clusters = sorted(taxon_df["cluster"].unique())
        genomes_comps = []
        for cluster_1, cluster_2 in itertools.combinations(sorted_clusters, 2):
            genomes_1 = taxon_df[taxon_df["cluster"] == cluster_1]["genome"]
            genomes_2 = taxon_df[taxon_df["cluster"] == cluster_2]["genome"]
            genomes_comps += list(itertools.product(genomes_1, genomes_2))
    else:
        species = taxon_df["clade"].unique()
        n_simu = taxon_df.shape[0] // len(species)
        genomes_comps = []
        for i in range(n_simu):
            genomes_comps += list(itertools.combinations([f"{spe}_{i}.fa" for spe in species], 2))
    return genomes_comps


def create_lastz_db(taxon_csv, genomes_path, cluster_name, db_name, threads, update=False, pairs=False, aligner="lastz"):
    """Creates or updates a sqlite3 database from a taxon csv file, generating all necessary alignments"""
    if not update:
        if os.path.exists(db_name):
            os.remove(db_name)
    sqlite3_conn = sqlite3.connect(db_name)
    taxon_df = pd.read_csv(taxon_csv)
    cur = sqlite3_conn.cursor()
    sorted_clusters = sorted(taxon_df[cluster_name].unique())
    if update:
        already_compared = cur.execute("SELECT genome1, genome2 FROM lastz").fetchall()
        # gather necessary genome comparisons
        genomes_comps = []
        for cluster_1, cluster_2 in itertools.combinations(sorted_clusters, 2):
            genomes_1 = sorted(taxon_df[taxon_df[cluster_name] == cluster_1]["genome"])
            genomes_2 = sorted(taxon_df[taxon_df[cluster_name] == cluster_2]["genome"])
            genomes_comps += list(itertools.product(genomes_1, genomes_2))

        sorted_genomes_comps = [(g1, g2) if g1 < g2 else (g2, g1) for g1, g2 in genomes_comps]
        genomes_comps = [genome_pair for genome_pair in sorted_genomes_comps if (genome_pair[0], genome_pair[1]) not in already_compared]
        genomes_comps = [(os.path.join(genomes_path, genome_1), os.path.join(genomes_path, genome_2)) for genome_1, genome_2 in genomes_comps]
        # print genome comps and already compared
        print(f"Already in database : {len(already_compared)} comparisons, {len(genomes_comps)} to align, total : {len(sorted_genomes_comps)} comparisons.")

    elif not update and not pairs:
        cur.execute("CREATE TABLE lastz (genome1 STRING, genome2 STRING, count_array blob, average_divergence REAL, total_matches INT, total_aligned INT, top_div REAL);")
        cur.execute("CREATE TABLE taxon (genome STRING, cluster STRING);")
        for _, row in taxon_df.iterrows():
            cur.execute("INSERT INTO taxon VALUES (?, ?)", (row["genome"], row[cluster_name]))

        # gather necessary genome comparisons
        genomes_comps = []
        for cluster_1, cluster_2 in itertools.combinations(sorted_clusters, 2):
            genomes_1 = taxon_df[taxon_df[cluster_name] == cluster_1]["genome"]
            genomes_2 = taxon_df[taxon_df[cluster_name] == cluster_2]["genome"]
            genomes_1 = [os.path.join(genomes_path, genome) for genome in sorted(genomes_1)]
            genomes_2 = [os.path.join(genomes_path, genome) for genome in sorted(genomes_2)]
            genomes_comps += list(itertools.product(genomes_1, genomes_2))

    elif not update and pairs:
        cur.execute("CREATE TABLE lastz (genome1 STRING, genome2 STRING, count_array blob, average_divergence REAL, total_matches INT, total_aligned INT, top_div REAL);")
        cur.execute("CREATE TABLE taxon (genome STRING, cluster STRING);")
        for _, row in taxon_df.iterrows():
            cur.execute("INSERT INTO taxon VALUES (?, ?)", (row["genome"], row[cluster_name]))
        genomes_comps = get_genome_comps(taxon_csv, pairs=True)
        genomes_comps = [[os.path.join(genomes_path, g1), os.path.join(genomes_path, g2)] for g1, g2 in genomes_comps]
    # run lastz in parallel
    batch_size = 2000
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        for i in range(0, len(genomes_comps), batch_size):
            # create unique prefix for mummer parallel runs
            num_genomes = len(genomes_comps[i:i+batch_size])
            prefixes = [f"mummer_{i}" for i in range(num_genomes)]
            res_list = list(executor.map(align_exec, genomes_comps[i:i+batch_size], [aligner]*num_genomes, prefixes))
            for res in res_list:
                align_entry(res, sqlite3_conn)
    return sqlite3_conn


def update_lastz_db(taxon_csv, genomes_path, cluster_name, db_name, threads):
    """Updates a sqlite3 database from a taxon csv file, generating all necessary alignments"""
    sqlite3_conn = sqlite3.connect(db_name)
    cur = sqlite3_conn.cursor()
    taxon_df = pd.read_csv(taxon_csv)
    sorted_clusters = sorted(taxon_df[cluster_name].unique())
    already_compared = cur.execute("SELECT genome1, genome2 FROM lastz").fetchall()
    # gather necessary genome comparisons
    genomes_comps = []
    for cluster_1, cluster_2 in itertools.combinations(sorted_clusters, 2):
        genomes_1 = sorted(taxon_df[taxon_df[cluster_name] == cluster_1]["genome"])
        genomes_2 = sorted(taxon_df[taxon_df[cluster_name] == cluster_2]["genome"])
        genomes_comps += list(itertools.product(genomes_1, genomes_2))
    sorted_genomes_comps = [(g1, g2) if g1 < g2 else (g2, g1) for g1, g2 in genomes_comps]
    genomes_comps = [genome_pair for genome_pair in sorted_genomes_comps if (genome_pair[0], genome_pair[1]) not in already_compared]
    genomes_comps = [(os.path.join(genomes_path, genome_1), os.path.join(genomes_path, genome_2)) for genome_1, genome_2 in genomes_comps]
    # print genome comps and already compared
    print(f"Already in database : {len(already_compared)} comparisons, {len(genomes_comps)} to align, total : {len(sorted_genomes_comps)} comparisons.")

    # run lastz in parallel
    # batching so that we can checkpoint
    batch_size = 2000
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        for i in range(0, len(genomes_comps), batch_size):
            print(f"Starting batch {i//batch_size + 1}/{len(genomes_comps)//batch_size + 1}")
            res_list = list(executor.map(align_exec, genomes_comps[i:i+batch_size]))
            for res in res_list:
                align_entry(res, sqlite3_conn)
    return sqlite3_conn
