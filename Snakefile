import pandas as pd


configfile: "config.yaml"

COHORT = config["cohort"]
SAMPLESHEET = pd.read_csv("samplesheet-" + COHORT + ".tsv", sep="\t", comment="#")
SAMPLESHEET = SAMPLESHEET.set_index('person', drop=False)

def metastatic_sites_option(person):
    if not pd.isna(SAMPLESHEET.loc[person, 'metastases']):
        sites = SAMPLESHEET.loc[person, 'metastases']
        return f'--metastatic-sites {sites}'
    else: 
        return ''

def full_citup_scan(samplesheet):
    samplesheet = samplesheet[~samplesheet.tests.isna()]
    short_scan = samplesheet[['person']].copy()
    short_scan['test'] = samplesheet.tests.apply(lambda x : [int(i) for i in x.split(',')])
    scan = short_scan.explode('test', ignore_index=True)
    return scan

def display_citup_scan(samplesheet):
    samplesheet = samplesheet[~samplesheet.tests.isna()]
    short_scan = samplesheet[['person']].copy()
    short_scan['display_tests'] = samplesheet.display_tests.apply(lambda x : [int(i) for i in x.split(',')])
    scan = short_scan.explode('display_tests', ignore_index=True)
    return scan

FULL_CITUP_SCAN = full_citup_scan(SAMPLESHEET)
DISPLAY_CITUP_SCAN = display_citup_scan(SAMPLESHEET)

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
        dagSolutions=lambda w: config["citup2"][w.cohort]["dag_solutions"],
        metastatic_sites = lambda w: metastatic_sites_option(w.person),
    resources:
        runtime=960,
        mem_mb=8 * 1024,
    log:
        "results/fixed_trees/{cohort}/{person}/logs/citup.log"
    threads: 8
    group:
        "citup2_fixed"
    shell:
        """
        echo "Job ID: ${{SLURM_JOB_ID:-NA}}" > {log}
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
            --mike-max-refined-solutions {params.dagSolutions} \\
            {params.args} --gurobi-threads {threads} >> {log} 2>&1
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
        dagSolutions=lambda w: config["citup2"][w.cohort]["dag_solutions"],
        metastatic_sites = lambda w: metastatic_sites_option(w.person),
    resources:
        runtime=48 * 60,
        mem_mb=8 * 1024,
    log:
        "results/penalty_trees/{cohort}/{person}/logs/citup.log"
    threads: 8
    group:
        "citup2_fixed"
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
            --numSolutions {params.numSolutions} \\
            --mike-max-refined-solutions {params.dagSolutions} \\
            {params.args} --gurobi-threads {threads} >> {log} 2>&1
        """


rule optimize_without_migrations_stage1:
    output:
        directory("results/optimal_trees/{cohort}/{person}/trees"),
    input:
        mutations="data/{cohort}/{person}/mutations.tsv",
    params:
        args=lambda w: config["citup2_no_migrations"][w.cohort]["args"],
        citup_py=lambda w: config["citup2_no_migrations"][w.cohort]["citup_py"],
        numSolutions=lambda w: config["citup2_no_migrations"][w.cohort]["num_solutions"],
        dagSolutions=lambda w: config["citup2_no_migrations"][w.cohort]["dag_solutions"],
    resources:
        runtime=24 * 60,
        mem_mb=8 * 1024,
    log:
        "results/optimal_trees/{cohort}/{person}/logs/citup.log"
    threads: 8
    group:
        "citup2_no_migations"
    shell:
        """
        echo "Job ID: ${{SLURM_JOB_ID:-NA}}" > {log}
        outdir={output[0]}
        mkdir -p $outdir
        outprefix="$outdir/{wildcards.person}"
        python {params.citup_py} \\
            -f {input.mutations} \\
            -o $outprefix \\
            --numSolutions {params.numSolutions} \\
            {params.args} --gurobi-threads {threads} >> {log} 2>&1
        """

rule migrations_for_citup_optima_stage1:
    output:
        directory("results/migration_trees/{cohort}/{person}/trees"),
    input:
        mutations="data/{cohort}/{person}/mutations.tsv",
        sites="data/{cohort}/{person}/sites.tsv",
        optima="results/optimal_trees/{cohort}/{person}/trees"
    params:
        args=lambda w: config["citup2"][w.cohort]["args"],
        citup_py=lambda w: config["citup2"][w.cohort]["citup_py"],
        numSolutions=lambda w: config["citup2"][w.cohort]["num_solutions"],
        dagSolutions=lambda w: config["citup2"][w.cohort]["dag_solutions"],
        metastatic_sites = lambda w: metastatic_sites_option(w.person),
        tree_scores = lambda w: expand("results/optimal_trees/{cohort}/{person}/trees/{person}.tree_0.txt",
                                        cohort=w.cohort, person=w.person),
    resources:
        runtime=960,
        mem_mb=8 * 1024,
    log:
        "results/migration_trees/{cohort}/{person}/logs/citup.log"
    threads: 8
    group:
        "citup2_no_migations"
    shell:
        """
        echo "Job ID: ${{SLURM_JOB_ID:-NA}}" > {log}
        outdir={output[0]}
        mkdir -p $outdir
        outprefix="$outdir/{wildcards.person}"

        python scripts/extract_fixed_tree.py {params.tree_scores} $outdir >> {log} 2>&1
        python {params.citup_py} \\
            -f {input.mutations} \\
            --sites-file {input.sites} \\
            --fixed-B-file $outdir/fixed_B.tsv \\
            --fixed-delta-file $outdir/fixed_delta.tsv {params.metastatic_sites} \\
            -o $outprefix \\
            --numSolutions {params.numSolutions} \\
            --mike-max-refined-solutions {params.dagSolutions} \\
            {params.args} --gurobi-threads {threads} >> {log} 2>&1
        """


rule min_migrations_for_citup_optima_stage1:
    output:
        directory("results/min_trees/{cohort}/{person}/trees"),
    input:
        mutations="data/{cohort}/{person}/mutations.tsv",
        sites="data/{cohort}/{person}/sites.tsv",
        #optima="results/optimal_trees/{cohort}/{person}/trees"
    params:
        args=lambda w: config["citup1c"][w.cohort]["args"],
        citup_py=lambda w: config["citup1c"][w.cohort]["citup_py"],
        numSolutions=lambda w: config["citup1c"][w.cohort]["num_solutions"],
        dagSolutions=lambda w: config["citup1c"][w.cohort]["dag_solutions"],
        metastatic_sites = lambda w: metastatic_sites_option(w.person),
        tree_scores = lambda w: expand("results/optimal_trees/{cohort}/{person}/trees/{person}.tree_0.txt",
                                        cohort=w.cohort, person=w.person),
    resources:
        runtime=990,
        mem_mb=8 * 1024,
    log:
        "results/min_trees/{cohort}/{person}/logs/citup.log"
    threads: 8
    shell:
        """
        echo "Job ID: ${{SLURM_JOB_ID:-NA}}" > {log}
        outdir={output[0]}
        mkdir -p $outdir
        outprefix="$outdir/{wildcards.person}"
        tree_score=$(python scripts/tree_score.py {params.tree_scores})
        python scripts/extract_fixed_tree.py {params.tree_scores} $outdir >> {log} 2>&1

        python {params.citup_py} \\
            -f {input.mutations} \\
            --sites-file {input.sites} \\
            -o $outprefix \\
            --fit-time-limit 57600 \\
            --warm-start-first-phase-b $outdir/fixed_B.tsv \\
            --warm-start-first-phase-delta $outdir/fixed_delta.tsv {params.metastatic_sites} \\
            --numSolutions {params.numSolutions} \\
            --minimize-migrations-at-score $tree_score --migration-score-epsilon 0.05 \\
            --migration-score-source without-reward \\
            {params.args} --gurobi-threads {threads} >> {log} 2>&1
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
        dagSolutions=lambda w: config["citup_scan"][w.cohort]["dag_solutions"],
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
            --mike-max-refined-solutions {params.dagSolutions} \\
            {params.args} --gurobi-threads {threads} >> {log} 2>&1
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
            --mike-max-refined-solutions {params.dagSolutions} \\
            --numSolutions {params.numSolutions} \\
            --mike-max-migrations {params.argmin} \\
            {params.args} --gurobi-threads {threads} >> {log} 2>&1
        """


rule run_conipher:
    input:
        "data/conipher/{case_id}/input.tsv"
    output:
        directory("results/conipher/{case_id}")
    resources:
        runtime=24 * 60,
        mem_mb=8 * 1024,
    shell:
        """
        Rscript CONIPHER-wrapper/src/run_treebuilding.R --input_tsv {input} \
            --min_cluster_size=3 --multi_trees=FALSE --pval_cutoff=0.05 \
            --out_dir {output} --prefix m --min_cluster_size=1
        """

rule conipher_trees_all:
    input:
        expand(
            "results/conipher/{person}",
            person=SAMPLESHEET.person,
        ),


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


rule stage1_optimal_trees_all:
    input:
        expand(
            "results/optimal_trees/{cohort}/{person}/trees",
            cohort=COHORT,
            person=SAMPLESHEET.person,
        ),


rule stage1_min_trees_all:
    input:
        expand(
            "results/min_trees/{cohort}/{person}/trees",
            cohort=COHORT,
            person=SAMPLESHEET.person,
        ),


rule stage1_migration_trees_all:
    input:
        expand(
            "results/migration_trees/{cohort}/{person}/trees",
            cohort=COHORT,
            person=SAMPLESHEET.person,
        ),


rule stage2_migration_scan_all:
    input:
        expand(
            "results/migration_scan/{cohort}/mig_{migrations:03d}/{person}",
            zip,
            cohort=[COHORT] * FULL_CITUP_SCAN.shape[0],
            person=FULL_CITUP_SCAN.person,
            migrations=FULL_CITUP_SCAN.test,
        ),


rule stage2_migration_scan_display_all:
    input:
        expand(
            "results/migration_scan/{cohort}/mig_{migrations:03d}/{person}",
            zip,
            cohort=[COHORT] * DISPLAY_CITUP_SCAN.shape[0],
            person=DISPLAY_CITUP_SCAN.person,
            migrations=DISPLAY_CITUP_SCAN.display_tests,
        ),


rule stage3_optimal_rerun_all:
    input:
        expand(
            "results/migration_scan/{cohort}/pareto/{person}/trees",
            cohort=COHORT,
            person=SAMPLESHEET.query("~argmin.isna()").person,
        ),
