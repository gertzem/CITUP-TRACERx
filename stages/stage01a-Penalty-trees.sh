set -euf

snakemake --profile snakemake_profile stage1_penalty_trees_all \
    --config cohort="$1"
