set -euf

echo "Running on the tests saved in samplesheet-$1.tsv"
snakemake --profile snakemake_profile stage2_migration_scan_all \
    --config cohort="$1"
