rule database:
    input:
        genomes=[os.path.join(config["genomes_dir"], g) for g in GENOME_LIST]
    params:
        up=lambda wc: "--update" if config["database_state"] == "update" else "",
        pairs=lambda wc: "--pairs" if config["pairs"]=="yes" else "",
        taxon_csv=config["taxon_csv"],
        genomes_dir=config["genomes_dir"],
        cluster_name=config["cluster_name"],
        aligner=config["aligner"]
    threads: config["max_threads"]
    output:
        database
    shell:
        "python lastz_parallel_db/main.py --threads {threads} --aligner {params.aligner} {params.up} {params.pairs} "
        "{params.taxon_csv} {params.genomes_dir} "
        "{params.cluster_name}  {output}"
