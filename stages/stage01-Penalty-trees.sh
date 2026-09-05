set -euf

snakemake --profile snakemake_profile stage1_penalty_trees_all \
    --config cohort="$1"

python scripts/read_citup_obj.py \ 
    $(find results/penalty_trees -name '*.tree_0.txt') > results/penalty_trees/tracerx/optima.tsv

python scripts/plot_by_tree_path.py \
    $(find results/penalty_trees -name trees)
