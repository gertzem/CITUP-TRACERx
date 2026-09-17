#! /bin/bash

set -euf


python scripts/read_citup_obj.py \
    $(find results/migration_trees -name '*.tree_0.txt') > "results/migration_trees/$1/optima.tsv"

python scripts/plot_by_tree_path.py \
    $(find "results/migration_trees/$1" -name trees)
