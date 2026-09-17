set -euf

snakemake --profile snakemake_profile stage1_migration_trees_all \
    --config cohort="$1"
