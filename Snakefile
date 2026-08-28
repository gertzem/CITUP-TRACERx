import pandas as pd


configfile: "config.yaml"


SAMPLESHEET = pd.read_csv("samplesheet.tsv", sep="\t", comment="#")


rule optimize_fixed_tree:
    output:
        directory("results/citup2_mets/{model}/{person}/trees"),
    input:
        mutations="data/tracerx/{person}/mutations.tsv",
        sites="data/tracerx/{person}/sites.tsv",
        B="data/tracerx/{person}/fixed_B.tsv",
        delta="data/tracerx/{person}/fixed_delta.tsv",
    params:
        args=lambda w: config["citup2"][w.model]["args"],
        citup_py=lambda w: config["citup2"][w.model]["citup_py"],
        numSolutions=lambda w: config["citup2"][w.model]["num_solutions"],
        metastatic_sites=lambda w: SAMPLESHEET.query(
            "person == @w.person"
        ).metastases.iloc[0],
    resources:
        runtime=960,
        mem_mb=8 * 1024,
    threads: 8
    group:
        "citup2_fixed"
    shell:
        """
        outdir={output[0]}
        mkdir -p $outdir
        outprefix="$outdir/{wildcards.person}"
        python {params.citup_py} \\
            -f {input.mutations} \\
            --sites-file {input.sites} \\
            --fixed-B-file {input.B} \\
            --fixed-delta-file {input.delta} \\
            --metastatic-sites {params.metastatic_sites} \\
            -o $outprefix \\
            --numSolutions {params.numSolutions} \\
            {params.args} --gurobi-threads {threads}
        """


rule rerun_with_optimal_migrations:
    output:
        directory("results/citup_scan_mets/{model}/pareto/{person}/trees"),
    input:
        mutations="data/tracerx/{person}/mutations.tsv",
        sites="data/tracerx/{person}/sites.tsv",
    params:
        args=lambda w: config["citup_optimal"][w.model]["args"],
        citup_py=lambda w: config["citup_optimal"][w.model]["citup_py"],
        numSolutions=lambda w: config["citup_optimal"][w.model]["num_solutions"],
        dagSolutions=lambda w: config["citup_optimal"][w.model]["dag_solutions"],
        metastatic_sites=lambda w: SAMPLESHEET.query(
            "person == @w.person"
        ).metastases.iloc[0],
        argmin=lambda w: int(SAMPLESHEET.query(
            "person == @w.person"
        ).argmin.iloc[0]),
    resources:
        runtime=24 * 60,
        mem_mb=8 * 1024,
    threads: 8
    shell:
        """
        outdir={output[0]}
        mkdir -p $outdir
        outprefix="$outdir/{wildcards.person}"
        python {params.citup_py} \\
            -f {input.mutations} \\
            --sites-file {input.sites} \\
            --metastatic-sites {params.metastatic_sites} \\
            -o $outprefix \\
            --enum-max-trees {params.dagSolutions} \\
            --numSolutions {params.numSolutions} \\
            --mike-max-migrations {params.argmin} \\
            {params.args} --gurobi-threads {threads}
        """



rule optimize_with_penalty:
    output:
        directory("results/citup2_free_mets/{model}/{person}/trees"),
    input:
        mutations="data/tracerx/{person}/mutations.tsv",
        sites="data/tracerx/{person}/sites.tsv",
    params:
        args=lambda w: config["citup2"][w.model]["args"],
        citup_py=lambda w: config["citup2"][w.model]["citup_py"],
        numSolutions=lambda w: config["citup2"][w.model]["num_solutions"],
        metastatic_sites=lambda w: SAMPLESHEET.query(
            "person == @w.person"
        ).metastases.iloc[0],
    resources:
        runtime=48 * 60,
        mem_mb=8 * 1024,
    threads: 8
    group:
        "citup2_fixed"
    shell:
        """
        outdir={output[0]}
        mkdir -p $outdir
        outprefix="$outdir/{wildcards.person}"
        python {params.citup_py} \\
            -f {input.mutations} \\
            --sites-file {input.sites} \\
            --metastatic-sites {params.metastatic_sites} \\
            -o $outprefix \\
            --numSolutions {params.numSolutions} \\
            {params.args} --gurobi-threads {threads}
        """


rule read_migration_scan:
    output:
        "results/citup_scan_mets/{method}/pareto/{person}/pareto.txt",
    input:
        lambda w: expand(
            "results/citup_scan_mets/{{method}}/mig_{migrations:03d}/{{person}}/{{person}}.mig_{migrations:03d}.treeScores.tsv",
            migrations=[int(m) for m in SAMPLESHEET.query("person==@w.person").tests.iloc[0].split(",")],
        ),
    localrule: True
    shell:
        """
        python scripts/pareto_front.py {input} > {output}
        """

rule optimize_for_fixed_migrations:
    output:
        "results/citup_scan_mets/{model}/mig_{migrations}/{person}/{person}.mig_{migrations}.treeScores.tsv",
    input:
        "data/tracerx/{person}/mutations.tsv",
        "data/tracerx/{person}/sites.tsv",
    params:
        args=lambda w: config["citup_scan"][w.model]["args"],
        citup_py=lambda w: config["citup_scan"][w.model]["citup_py"],
        numSolutions=lambda w: config["citup_scan"][w.model]["num_solutions"],
        metastatic_sites=lambda w: SAMPLESHEET.query(
            "person == @w.person"
        ).metastases.iloc[0],
        max_migrations=lambda w: int(w.migrations),
    resources:
        runtime=24 * 60,
        mem_mb=8 * 1024,
    threads: 8
    shell:
        """
        output={output}
        prefix=${{output%.treeScores.tsv}}
        python {params.citup_py} \\
            -f {input[0]} \\
            --sites-file {input[1]} \\
            -o $prefix \\
            --metastatic-sites {params.metastatic_sites} \\
            --numSolutions {params.numSolutions} \\
            --mike-max-migrations {params.max_migrations} \\
            {params.args} --gurobi-threads {threads}
        """

rule plot_for_fixed_migrations:
    output:
        directory("results/citup_scan_mets/{model}/mig_{migrations}/{person}/plots"),
    input:
        "results/citup_scan_mets/{model}/mig_{migrations}/{person}/{person}.mig_{migrations}.treeScores.tsv",
        "data/tracerx/{person}/sites.tsv",
    localrule: True
    shell:
        """
        mkdir -p {output}
        indir=$(dirname {input[0]})

        python ../CITUP/plot_with_gene_names.py "$indir" \\
            --output-dir {output} \\
            --sample-sites {input[1]} \\
            --with-site-transitions \\
            --with_node_info \\
            --all-mike-dag-colorings
        """

rule plot_for_fixed_tree:
    output:
        directory("results/citup2_mets/{model}/{person}/plots"),
    input:
        "results/citup2_mets/{model}/{person}/trees",
        "data/tracerx/{person}/sites.tsv",
    localrule: True
    shell:
        """
        mkdir -p {output}
        indir={input[0]}

        python ../CITUP/plot_with_gene_names.py "$indir" \\
            --output-dir {output} \\
            --sample-sites {input[1]} \\
            --with-site-transitions \\
            --with_node_info \\
            --all-mike-dag-colorings
        """


rule read_fixed_tree_optima:
    output: 
        "results/citup2_mets/conipher/optima.tsv",
    input:
        expand("results/citup2_mets/conipher/{person}/trees/{person}.tree_0.txt", person=SAMPLESHEET.person),
    shell:
        """
        python scripts/read_citup_obj.py {input} > {output}
        """


rule read_optima_with_penalty:
    output: 
        "results/citup2_free_mets/conipher/optima.tsv",
    input:
        expand("results/citup2_free_mets/conipher/{person}/trees/{person}.tree_0.txt", 
                person=SAMPLESHEET.query('~free_mets_optima.isna()').person),
    shell:
        """
        python scripts/read_citup_obj.py {input} > {output}
        """


rule perform_migration_scan:
    input:
        expand("results/citup_scan_mets/conipher/pareto/{person}/trees", 
                person=SAMPLESHEET.query('~argmin.isna()').person),


rule rescore_fixed_trees:
    input:
        expand("results/citup2_mets/conipher/{person}/trees", person=SAMPLESHEET.person),



rule build_penalty_trees:
    input:
        expand(
            "results/citup2_free_mets/conipher/{person}/trees",
            person=SAMPLESHEET.person,
        ),


rule citup_conipher_pareto_mets:
    input:
        expand(
            "results/citup_scan_mets/conipher/pareto/{person}/pareto.txt",
            person=SAMPLESHEET.person,
        ),


rule read_pareto_fronts:
    output:
        "results/citup_scan_mets/conipher/pareto/optima.tsv"
    input: 
        expand(
            "results/citup_scan_mets/conipher/pareto/{person}/pareto.txt",
            person=SAMPLESHEET.person,
        ),
    params:
        samples=lambda w: SAMPLESHEET.to_json(orient="records"),
    script:
        "read_pareto_optimal.py"
