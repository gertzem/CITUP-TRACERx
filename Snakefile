import pandas as pd


configfile: "config.yaml"


SAMPLESHEET = pd.read_csv(config["samplesheet"], sep="\t", comment="#")
SAMPLESHEET = SAMPLESHEET.set_index('person', drop=False)
COHORT = config["cohort"]

def metastatic_sites_option(person):
    if not pd.isna(SAMPLESHEET.loc[person, 'metastases']):
        sites = SAMPLESHEET.loc[person, 'metastases']
        return f'--metastatic-sites {sites}'
    else: 
        return ''

rule optimize_fixed_tree_stage0:
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


rule optimize_with_penalty_stage1:
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


rule scan_for_target_migrations_stage2:
    output:
        directory("results/migration_scan/{cohort}/mig_{migrations}/{person}")
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
    log:
        "results/migration_scan/{cohort}/logs/mig_{migrations}/{person}/citup.log"
    threads: 8
    shell:
        """
        echo "Job ID: ${{SLURM_JOB_ID:-NA}}" > {log}
        mkdir -p {output}
        prefix={output}/{wildcards.person}.mig_{wildcards.migrations}
        python {params.citup_py} \\
            -f {input[0]} \\
            --sites-file {input[1]} \\
            -o $prefix {params.metastatic_sites} \\
            --numSolutions {params.numSolutions} \\
            --mike-max-migrations {params.max_migrations} \\
            {params.args} --gurobi-threads {threads} >> {log} 2>&1
        """


rule read_pareto_front_stage2:
    output:
        "results/migration_scan/{method}/pareto/{person}/pareto.txt",
    input:
        lambda w: expand(
            "results/migration_scan/{{method}}/mig_{migrations:03d}/{{person}}",
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


rule rerun_with_optimal_migrations_stage3:
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
    log:
        "results/migration_scan/{cohort}/pareto/logs/{person}/citup.log"
    shell:
        """
        echo "Job ID: ${{SLURM_JOB_ID:-NA}}" > {log}
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
            {params.args} --gurobi-threads {threads} >> {log} 2>&1
        """



rule stage0_fixed_trees_all:
    input:
        expand(
            "results/fixed_trees/{cohort}/{person}/trees",
            cohort=COHORT,
            person=SAMPLESHEET.person,
        ),


rule stage1_penalty_trees_all:
    input:
        expand(
            "results/penalty_trees/{cohort}/{person}/trees",
            cohort=COHORT,
            person=SAMPLESHEET.person,
        ),


rule stage2_migration_scan_all:
    input:
        expand(
            "results/migration_scan/{cohort}/pareto/{person}/pareto.txt",
            cohort=COHORT,
            person=SAMPLESHEET.person,
        ),


rule stage3_optimal_rerun_all:
    input:
        expand(
            "results/migration_scan/{cohort}/pareto/{person}/trees",
            cohort=COHORT,
            person=SAMPLESHEET.query("~argmin.isna()").person,
        ),
