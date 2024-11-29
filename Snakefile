configfile: "./smk_config.yml"
import itertools
import os

import pandas as pd


CLUSTER_LIST = sorted(list(set(pd.read_csv(config["taxon_csv"])[config["cluster_name"]])))
GENOME_LIST = sorted(list(set(pd.read_csv(config["taxon_csv"])["genome"])))

# basic rules
database = os.path.join(config["results_dir"], config["database_name"])
plot=[f"{config['results_dir']}fig2_plots/{bac1}_vs_{bac2}_plot_fig2.png" for bac1, bac2 in itertools.combinations(CLUSTER_LIST, 2)]
full_mld=[f"{config['results_dir']}full_mlds/{bac1}_vs_{bac2}_full_mld_comp.csv" for bac1, bac2 in itertools.combinations(CLUSTER_LIST, 2)]
fitted_params=[f"{config['results_dir']}fitted_params/{bac1}_vs_{bac2}_fitted_params.csv" for bac1, bac2 in itertools.combinations(CLUSTER_LIST, 2)]
binned_mld=[f"{config['results_dir']}binned_mlds/{bac1}_vs_{bac2}_binned_mld.csv" for bac1, bac2 in itertools.combinations(CLUSTER_LIST, 2)]
lengths_a=[config['results_dir'] + len_distr for len_distr in expand("lengths_distributions/{fasta_dir}_distribution.{ext}", fasta_dir = CLUSTER_LIST, ext = ["png", "csv"])]
surfaces=[f"{config['results_dir']}surfaces/{bac1}_vs_{bac2}_surface_plot.png" for bac1, bac2 in itertools.combinations(CLUSTER_LIST, 2)]
L0s=f"{config['results_dir']}all_L0s.csv"
tree=config["results_dir"] + config["tree_annotation"] + "_tree_big.svg"

rule_all_list = [plot, full_mld, fitted_params, binned_mld, lengths_a, surfaces, L0s, tree]

# genome wise fits
comparisons=[f"{config['results_dir']}analyse_comparisons/{bac1}_vs_{bac2}_inflexion_res.csv" for bac1, bac2 in itertools.combinations(CLUSTER_LIST, 2)]



if config["from_f_mld"] == "yes":
    include: "florian_mld.smk"
elif config["database_state"] != "complete":
    rule_all_list.append(database)
    include: "database.smk"

if config["genome_wise_fit"] == "yes":
    rule_all_list.extend(comparisons)
    include: "genome_wise_inflexion.smk"
else:
    include: "cluster_wise_inflexion.smk"

rule all:
    input:
        rule_all_list


rule pretreat:
    output:
        genomes=temp([os.path.join(config["pretreatment_dir"], g) for g in GENOME_LIST])
    params:
        genome_dir=config["genomes_dir"],
        pretreatment_dir=config["pretreatment_dir"],
        lastz_tools_dir=config["lastz_tools_dir"]
    threads: config["max_threads"]
    shell:
        "python pretreatment/pretreat.py --threads {threads} {params.lastz_tools_dir} {params.genome_dir} {params.pretreatment_dir}"


rule merge:
    input:
        database_path=database
    output:
        full_mld=[f"{config['results_dir']}full_mlds/{bac1}_vs_{bac2}_full_mld_comp.csv" for bac1, bac2 in itertools.combinations(CLUSTER_LIST, 2)],
        binned_mld=[f"{config['results_dir']}binned_mlds/{bac1}_vs_{bac2}_binned_mld.csv" for bac1, bac2 in itertools.combinations(CLUSTER_LIST, 2)]
    params:
        taxon_csv=config["taxon_csv"],
        cluster_name=config["cluster_name"],
        full_mld_dir=config["results_dir"] + "full_mlds/",
        binned_mld_dir=config["results_dir"] + "binned_mlds/",
        pairs=lambda wc: "--pairs" if config["pairs"]=="yes" else ""

    threads: config["max_threads"]
    shell:
        "python parse/main.py {params.pairs} --threads {threads} --from_sqlite_db {input.database_path} --full_mld {params.full_mld_dir} "
        "--binned_mld {params.binned_mld_dir} --taxon_csv {params.taxon_csv} --cluster_name {params.cluster_name}"


rule plot:
    input:
        binned_mld=config["results_dir"] + "binned_mlds/{species_1}_vs_{species_2}_binned_mld.csv",
        fitted_params=config["results_dir"] + "fitted_params/{species_1}_vs_{species_2}_fitted_params.csv"
    output:
        config["results_dir"] + "fig2_plots/{species_1}_vs_{species_2}_plot_fig2.png"
    params:
        species=lambda w: f"{w.species_1},{w.species_2}",
        results_dir=config["results_dir"],
        mus=config["mus"],
        muc=config["muc"],
        delta=config["delta"]
    script:
        "plot_mld_fit.R"


rule lengths:
    input:
        config["genomes_dir"]
    output:
        csv_distr=expand(config["results_dir"] + "lengths_distributions/{cluster}_distribution.csv", cluster = CLUSTER_LIST),
        histo=expand(config["results_dir"] + "lengths_distributions/{cluster}_distribution.png", cluster = CLUSTER_LIST)
    params:
        taxon_csv=config["taxon_csv"],
        cluster_label=config["cluster_name"],
        output_dir=config["results_dir"] + "lengths_distributions/"
    shell:
        "python length_analysis.py --taxon_csv {params.taxon_csv} --cluster_label {params.cluster_label} --save_dir {params.output_dir} {input}"


rule L0:
    input:
        distribs=[config['results_dir'] + len_distr for len_distr in expand("lengths_distributions/{cluster}_distribution.csv", cluster = CLUSTER_LIST)]
    output:
        config["results_dir"] + "all_L0s.csv"
    params:
        distr_dir=config["results_dir"] + "lengths_distributions/"
    shell:
        "python get_L0.py {params.distr_dir} {output}"


rule fit:
    input:
        binned_mld=config["results_dir"] + "binned_mlds/{species_1}_vs_{species_2}_binned_mld.csv",
        all_L0s=config["results_dir"] + "all_L0s.csv"
    output:
        fitted_params=config["results_dir"] + "fitted_params/{species_1}_vs_{species_2}_fitted_params.csv",
        surface_plot=config["results_dir"] + "surfaces/{species_1}_vs_{species_2}_surface_plot.png"
    params:
        species=lambda w: f"{w.species_1},{w.species_2}",
        mus=config["mus"],
        muc=config["muc"],
        delta=config["delta"]
    shell:
        "python fit/main.py --bacs '{params.species}'  --L0 {input.all_L0s} "
        "--mus {params.mus} --muc {params.muc} --delta {params.delta} "
        "--save_surface_plot '{output.surface_plot}' '{input.binned_mld}' '{output.fitted_params}'"


rule trees:
    input:
        fitted_params=[f"{config['results_dir']}fitted_params/{bac1}_vs_{bac2}_fitted_params.csv" for bac1, bac2 in itertools.combinations(CLUSTER_LIST, 2)],
        lengths=[f"{config['results_dir']}lengths_distributions/{bac}_distribution.csv" for bac in CLUSTER_LIST],
        inflexion_file=config["results_dir"] + "inflexion_exists.csv",
        inflexion_percentage=config['results_dir'] + "inflexion_by_cluster.csv"
    output:
        config["results_dir"] + config["tree_annotation"] + "_tree_big.svg",
        config["results_dir"] + "fitteddistance_vs_founddistance.png",
        config["results_dir"] + "hist_fitteddistance.png"
    params:
        taxon_csv=config["taxon_csv"],
        cluster_name=config["cluster_name"],
        tree_annotation=config["tree_annotation"],
        fitted_params_dir=config["results_dir"] + "fitted_params/",
        results_dir=config["results_dir"],
        genome_wise_fit=config["genome_wise_fit"],
        genome_lengths=config["results_dir"] + "lengths_distributions/",
        filter_min_genomes=config["filter_min_genomes"]
    script:
        "make_trees.R"
