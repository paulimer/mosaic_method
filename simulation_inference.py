#!/usr/bin/env python3

# imports
import os

import yaml

from parse.fun import *
from lastz_parallel_db.utils import *
from fit.fit import *
from plot_mld_fit import *

configfile = "inf_sim_cfg.yml"
with open(configfile, 'r') as stream:
    cfg = yaml.safe_load(stream)
os.makedirs(cfg["results_dir"], exist_ok=True)

# running inference ------------------------------------------------------------
# alignment
database_path = os.path.join(cfg["results_dir"], cfg["database_name"])
con = create_lastz_db(
    cfg["taxon_csv"],
    cfg["genomes_dir"],
    "clade",
    database_path,
    cfg["max_threads"],
    False,
    False,
    cfg["aligner"]
)
con.close()

# mlds
taxon_df = pd.read_csv(cfg["taxon_csv"])
level_list = sorted(taxon_df["clade"].unique())
levels = list(itertools.combinations(level_list, 2))
binned_mlds = {}
for level in levels:
    genome_comps = get_genome_comp(level, cfg["taxon_csv"], "", "clade", output_csv=False, pairs=False)
    full_mld = get_all_mlds(genome_comps, database_path, threads=cfg["max_threads"])
    summed_mld = sum_mlds(full_mld)
    binned_mld = bin_mld(
        summed_mld,
        linear_bin_width=3,
        limit_size=30.5,
        power_increment=0.1,
        ncomp=full_mld.shape[0]
    )
    binned_mlds[level] = binned_mld

# fits
res_opt_full_dict = {}
res_opt_minus3_dict = {}
for level in levels:
    res_opt_full, res_opt_minus3 = fit_params(
        "dual-annealing",
        np.array([5, -5]),
        binned_mlds[level]["freq"],
        0.1,
        binned_mlds[level]["match_length"],
        cfg["mus"],
        cfg["muc"],
        cfg["delta"],
        float(cfg["megagene_size"])*cfg["n_megagenes"],
    )
    res_opt_full_dict[level] = res_opt_full
    res_opt_minus3_dict[level] = res_opt_minus3


# plot
os.makedirs(os.path.join(cfg["results_dir"], "plots"), exist_ok=True)

for level in levels:
    plot_mld_fit(
        binned_mlds[level],
        cfg["muc"],
        cfg["mus"],
        cfg["delta"],
        res_opt_full_dict[level].x,
        level,
        float(cfg["megagene_size"])*cfg["n_megagenes"],
        os.path.join(cfg["results_dir"], "plots", f"{level[0]}_{level[1]}.png")
    )


# saving output ----------------------------------------------------------------
os.makedirs(os.path.join(cfg["results_dir"], "binned_mlds"), exist_ok=True)

for level in levels:
    binned_mlds[level].to_csv(os.path.join(cfg["results_dir"], "binned_mlds", f"{level[0]}_{level[1]}.csv"), index=False)

res_list = []
for level in levels:
    res_list.append({
        "level": f"{level[0]}_{level[1]}",
        "log10tau": res_opt_full_dict[level].x[0],
        "log10rho": res_opt_full_dict[level].x[1],
        "minimum": res_opt_full_dict[level].fun,
        "minimum_minus3": res_opt_minus3_dict[level].fun
    })
res_df = pd.DataFrame(res_list)
res_df.to_csv(os.path.join(cfg["results_dir"], "fit_results.csv"), index=False)
