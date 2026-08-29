import pandas as pd


configfile: "config.yaml"


SAMPLESHEET = pd.read_csv(config["samplesheet"], sep="\t", comment="#")
SAMPLESHEET = SAMPLESHEET.set_index('person', drop=False)
COHORT = config["cohort"]

def metastatic_sites_option(person):
    sites = SAMPLESHEET.loc[person, 'metastases']
    if sites:
        return f'--metastatic-sites {sites}'
    else: 
        return ''

rule optimize_fixed_tree:
    output:
        directory("results/fixed_trees/{cohort}/{person}/trees"),
    input:
        mutations="data/{cohort}/{person}/mutations.tsv",
        sites="data/{cohort}/{person}/sites.tsv",
        B="data/{cohort}/{person}/fixed_B.tsv",
        delta="data/{cohort}/{person}/fixed_delta.tsv",
    params:
        args=lambda w: config["citup2"][w.cohort]["args"],
        citup_py=lambda w: config["citup2"][w.cohort]["citup_py"],
        numSolutions=lambda w: config["citup2"][w.cohort]["num_solutions"],
        metastatic_sites = lambda w: metastatic_sites_option(w.person),
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
            --fixed-delta-file {input.delta} {params.metastatic_sites} \\
            -o $outprefix \\
            --numSolutions {params.numSolutions} \\
            {params.args} --gurobi-threads {threads}
        """


rule read_fixed_tree_optima:
    output:
        "results/fixed_trees/{cohort}/optima.tsv",
    input:
        lambda w: expand(
            "results/fixed_trees/{w.cohort}/{person}/trees/{person}.tree_0.txt",
            person=SAMPLESHEET.person,
        ),
    shell:
        """
        python scripts/read_citup_obj.py {input} > {output}
        """


rule plot_for_fixed_tree:
    output:
        directory("results/fixed_trees/{cohort}/{person}/plots"),
    input:
        "results/fixed_trees/{cohort}/{person}/trees",
        "data/{cohort}/{person}/sites.tsv",
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


rule optimize_with_penalty:
    output:
        directory("results/penalty_trees/{cohort}/{person}/trees"),
    input:
        mutations="data/{cohort}/{person}/mutations.tsv",
        sites="data/{cohort}/{person}/sites.tsv",
    params:
        args=lambda w: config["citup2"][w.cohort]["args"],
        citup_py=lambda w: config["citup2"][w.cohort]["citup_py"],
        numSolutions=lambda w: config["citup2"][w.cohort]["num_solutions"],
        metastatic_sites = lambda w: metastatic_sites_option(w.person),
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
            --sites-file {input.sites} {params.metastatic_sites} \\
            -o $outprefix \\
            --numSolutions {params.numSolutions} \\
            {params.args} --gurobi-threads {threads}
        """


rule read_optima_with_penalty:
    output:
        "results/penalty_trees/{cohort}/optima.tsv",
    input:
        lambda w:
        expand(
            "results/penalty_trees/{w.cohort}/{person}/trees/{person}.tree_0.txt",
            person=SAMPLESHEET.query("~free_mets_optima.isna()").person,
        ),
    shell:
        """
        python scripts/read_citup_obj.py {input} > {output}
        """


rule scan_for_target_migrations:
    output:
        "results/migration_scan/{cohort}/mig_{migrations}/{person}/{person}.mig_{migrations}.treeScores.tsv",
    input:
        "data/{cohort}/{person}/mutations.tsv",
        "data/{cohort}/{person}/sites.tsv",
    params:
        args=lambda w: config["citup_scan"][w.cohort]["args"],
        citup_py=lambda w: config["citup_scan"][w.cohort]["citup_py"],
        numSolutions=lambda w: config["citup_scan"][w.cohort]["num_solutions"],
        metastatic_sites = lambda w: metastatic_sites_option(w.person),
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
            -o $prefix {params.metastatic_sites} \\
            --numSolutions {params.numSolutions} \\
            --mike-max-migrations {params.max_migrations} \\
            {params.args} --gurobi-threads {threads}
        """


rule plot_for_target_migrations:
    output:
        directory("results/migration_scan/{cohort}/mig_{migrations}/{person}/plots"),
    input:
        "results/migration_scan/{cohort}/mig_{migrations}/{person}/{person}.mig_{migrations}.treeScores.tsv",
        "data/{cohort}/{person}/sites.tsv",
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


rule read_migration_scan:
    output:
        "results/migration_scan/{method}/pareto/{person}/pareto.txt",
    input:
        lambda w: expand(
            "results/migration_scan/{{method}}/mig_{migrations:03d}/{{person}}/{{person}}.mig_{migrations:03d}.treeScores.tsv",
            migrations=[
                int(m)
                for m in SAMPLESHEET.query("person==@w.person")
                .tests.iloc[0]
                .split(",")
            ],
        ),
    localrule: True
    shell:
        """
        python scripts/pareto_front.py {input} > {output}
        """


rule rerun_with_optimal_migrations:
    output:
        directory("results/migration_scan/{cohort}/pareto/{person}/trees"),
    input:
        mutations="data/{cohort}/{person}/mutations.tsv",
        sites="data/{cohort}/{person}/sites.tsv",
    params:
        args=lambda w: config["citup_optimal"][w.cohort]["args"],
        citup_py=lambda w: config["citup_optimal"][w.cohort]["citup_py"],
        numSolutions=lambda w: config["citup_optimal"][w.cohort]["num_solutions"],
        dagSolutions=lambda w: config["citup_optimal"][w.cohort]["dag_solutions"],
        metastatic_sites = lambda w: metastatic_sites_option(w.person),
        argmin=lambda w: int(SAMPLESHEET.query("person == @w.person").argmin.iloc[0]),
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
            --sites-file {input.sites} {params.metastatic_sites} \\
            -o $outprefix \\
            --enum-max-trees {params.dagSolutions} \\
            --numSolutions {params.numSolutions} \\
            --mike-max-migrations {params.argmin} \\
            {params.args} --gurobi-threads {threads}
        """


rule optimal_rerun_all:
    input:
        expand(
            "results/migration_scan/{cohort}/pareto/{person}/trees",
            cohort=COHORT,
            person=SAMPLESHEET.query("~argmin.isna()").person,
        ),


rule fixed_trees_all:
    input:
        expand(
            "results/fixed_trees/{cohort}/{person}/trees",
            cohort=COHORT,
            person=SAMPLESHEET.person,
        ),


rule penalty_trees_all:
    input:
        expand(
            "results/penalty_trees/{cohort}/{person}/trees",
            cohort=COHORT,
            person=SAMPLESHEET.person,
        ),


rule pareto_fronts_all:
    input:
        expand(
            "results/migration_scan/{cohort}/pareto/{person}/pareto.txt",
            cohort=COHORT,
            person=SAMPLESHEET.person,
        ),


rule read_pareto_optima_all:
    output:
        "results/migration_scan/{cohort}/pareto/optima.tsv",
    input:
        expand(
            "results/migration_scan/{cohort}/pareto/{person}/pareto.txt",
            cohort=COHORT,
            person=SAMPLESHEET.person,
        ),
    params:
        samples=lambda w: SAMPLESHEET.to_json(orient="records"),
    script:
        "read_pareto_optimal.py"
