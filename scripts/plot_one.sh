#! /bin/bash

set -euf

indir=$1
person=$(basename $indir)
outdir=$(dirname $(dirname $indir))/plots

python ../CITUP/plot_with_gene_names.py "$indir" \
            --output-dir "$outdir" \
            --sample-sites "data/$cohort/$person/sites.tsv" \
            --with-site-transitions \
            --with_node_info \
            --all-mike-dag-colorings

python ../CITUP/plot_commigration_graphs.py "$indir" \
    --output-dir "$outdir" \
    --sites-file "data/$cohort/$person/sites.tsv"