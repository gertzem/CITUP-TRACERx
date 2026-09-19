set -euf

snakemake --profile snakemake_profile stage0_fixed_trees_all \
    --config cohort="$1"

python scripts/read_citup_obj.py $(find "results/fixed_trees/$1" -name '*.tree_0.txt') | sort \
    > "results/fixed_trees/$1/optima.tsv"

python scripts/plot_by_tree_path.py \
    $(find "results/fixed_trees/$1" -name trees)
