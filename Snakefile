import pandas as pd

configfile: "config.yaml"

SAMPLESHEET = pd.read_csv('samplesheet.tsv', sep="\t")

rule citup2:
    output:
        directory("results/citup2/{model}/{person}/trees"),
    input:
        "data/tracerx/{person}/mutations.tsv",
        "data/tracerx/{person}/sites.tsv"
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
        mutations = "data/tracerx/{person}/mutations.tsv",
        sites = "data/tracerx/{person}/sites.tsv",
        B = "data/tracerx/{person}/fixed_B.tsv",
        delta = "data/tracerx/{person}/fixed_delta.tsv"
    params:
        args=lambda w: config["citup2"][w.model]["args"],
        citup_py=lambda w: config["citup2"][w.model]["citup_py"],
        numSolutions=lambda w: config["citup2"][w.model]["num_solutions"],
    resources:
        runtime=960,
        mem_mb=8 * 1024,
    threads: 8
    group: "citup2_fixed",
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
        mutations = "data/tracerx/{person}/mutations.tsv",
        sites = "data/tracerx/{person}/sites.tsv",
        B = "data/tracerx/{person}/fixed_B.tsv",
        delta = "data/tracerx/{person}/fixed_delta.tsv"
    params:
        args=lambda w: config["citup2"][w.model]["args"],
        citup_py=lambda w: config["citup2"][w.model]["citup_py"],
        numSolutions=lambda w: config["citup2"][w.model]["num_solutions"],
        metastatic_sites=lambda w: SAMPLESHEET.query('person == @w.person').metastases.iloc[0],
    resources:
        runtime=960,
        mem_mb=8 * 1024,
    threads: 8
    group: "citup2_fixed",
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

rule citup2_free_mets:
    output:
        directory("results/citup2_free_mets/{model}/{person}/trees"),
    input:
        mutations = "data/tracerx/{person}/mutations.tsv",
        sites = "data/tracerx/{person}/sites.tsv",
    params:
        args=lambda w: config["citup2"][w.model]["args"],
        citup_py=lambda w: config["citup2"][w.model]["citup_py"],
        numSolutions=lambda w: config["citup2"][w.model]["num_solutions"],
        metastatic_sites=lambda w: SAMPLESHEET.query('person == @w.person').metastases.iloc[0],
    resources:
        runtime=960,
        mem_mb=8 * 1024,
    threads: 8
    group: "citup2_fixed",
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


rule citup_conipher_free:
    input:
        expand("results/citup2/conipher/{person}/trees", person=SAMPLESHEET.person)

rule citup_conipher_fixed:
    input:
        expand("results/citup2_fixed/conipher/{person}/trees", person=SAMPLESHEET.person)


rule citup_conipher_mets:
    input:
        expand("results/citup2_mets/conipher/{person}/trees", person=SAMPLESHEET.person)


rule citup_conipher_free_mets:
    input:
        expand("results/citup2_free_mets/conipher/{person}/trees", person=SAMPLESHEET.person)

