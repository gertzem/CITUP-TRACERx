import pandas as pd


configfile: "config.yaml"


SAMPLESHEET = pd.read_csv("samplesheet.tsv", sep="\t", comment="#")


rule citup2:
    output:
        directory("results/citup2/{model}/{person}/trees"),
    input:
        "data/tracerx/{person}/mutations.tsv",
        "data/tracerx/{person}/sites.tsv",
    params:
        args=lambda w: config["citup2"][w.model]["args"],
        citup_py=lambda w: config["citup2"][w.model]["citup_py"],
        numSolutions=lambda w: config["citup2"][w.model]["num_solutions"],
    resources:
        runtime=960,
        mem_mb=8 * 1024,
    threads: 8
    shell:
        """
        outdir={output[0]}
        mkdir -p $outdir
        outprefix="$outdir/{wildcards.person}"
        python {params.citup_py} \
            -f {input[0]} \
            --sites-file {input[1]} \
            -o $outprefix \
            --numSolutions {params.numSolutions} \
            {params.args} --gurobi-threads {threads}
        """


rule citup2_fixed:
    output:
        directory("results/citup2_fixed/{model}/{person}/trees"),
    input:
        mutations="data/tracerx/{person}/mutations.tsv",
        sites="data/tracerx/{person}/sites.tsv",
        B="data/tracerx/{person}/fixed_B.tsv",
        delta="data/tracerx/{person}/fixed_delta.tsv",
    params:
        args=lambda w: config["citup2"][w.model]["args"],
        citup_py=lambda w: config["citup2"][w.model]["citup_py"],
        numSolutions=lambda w: config["citup2"][w.model]["num_solutions"],
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
            -o $outprefix \\
            --numSolutions {params.numSolutions} \\
            {params.args} --gurobi-threads {threads}
        """


rule citup2_mets:
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


rule citup2_optimal_mets:
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
        runtime=960,
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



rule citup2_free_mets:
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
            --metastatic-sites {params.metastatic_sites} \\
            -o $outprefix \\
            --numSolutions {params.numSolutions} \\
            {params.args} --gurobi-threads {threads}
        """


rule citup_pareto:
    output:
        "results/citup_scan/{method}/{person}/trees/pareto.txt",
    input:
        lambda w: expand(
            "results/citup_scan/{{method}}/{{person}}/trees/mig{migrations}.treeScores.tsv",
            migrations=SAMPLESHEET.query("person==@w.person").tests.iloc[0].split(","),
        ),
    localrule: True
    shell:
        """
        python scripts/pareto_front.py {input} > {output}
        """


rule citup_migration:
    output:
        "results/citup_scan/{model}/{person}/trees/mig{migrations}.treeScores.tsv",
    input:
        "data/tracerx/{person}/mutations.tsv",
        "data/tracerx/{person}/sites.tsv",
    params:
        args=lambda w: config["citup_scan"][w.model]["args"],
        citup_py=lambda w: config["citup_scan"][w.model]["citup_py"],
        numSolutions=lambda w: config["citup_scan"][w.model]["num_solutions"],
    resources:
        runtime=960,
        mem_mb=16 * 1024,
    threads: 8
    shell:
        """
        output={output}
        prefix=${{output%.treeScores.tsv}}
        python {params.citup_py} \
            -f {input[0]} \
            --sites-file {input[1]} \
            -o $prefix \
            --numSolutions {params.numSolutions} \
            --mike-max-migrations {wildcards.migrations} \
            {params.args} --gurobi-threads {threads}
        """

rule citup_pareto_mets:
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

rule citup_migration_mets:
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

rule citup_conipher_free:
    input:
        expand("results/citup2/conipher/{person}/trees", person=SAMPLESHEET.person),


rule citup_mets_optimal:
    input:
        expand("results/citup_scan_mets/conipher/pareto/{person}/trees", 
                person=SAMPLESHEET.query('~argmin.isna()').person),


rule citup_conipher_fixed:
    input:
        expand(
            "results/citup2_fixed/conipher/{person}/trees", person=SAMPLESHEET.person
        ),


rule citup_conipher_mets:
    input:
        expand("results/citup2_mets/conipher/{person}/trees", person=SAMPLESHEET.person),


rule citup_conipher_free_mets:
    input:
        expand(
            "results/citup2_free_mets/conipher/{person}/trees",
            person=SAMPLESHEET.person,
        ),


rule citup_conipher_pareto:
    input:
        expand(
            "results/citup_scan/conipher/{person}/trees/pareto.txt",
            person=SAMPLESHEET.person,
        ),

rule citup_mig8a_pareto:
    input:
        expand(
            "results/citup_scan/mig8a/{person}/trees/pareto.txt",
            person=SAMPLESHEET.person,
        ),


rule citup_mig8_pareto:
    input:
        expand(
            "results/citup_scan/mig8/{person}/trees/pareto.txt",
            person=SAMPLESHEET.person,
        ),


rule citup_conipher_pareto_mets:
    input:
        expand(
            "results/citup_scan_mets/conipher/pareto/{person}/pareto.txt",
            person=SAMPLESHEET.person,
        ),
