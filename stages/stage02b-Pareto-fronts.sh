set -euf

bash scripts/all-pareto-fronts.sh "samplesheet-$1.tsv" "results/migration_scan/$1"
python scripts/pareto_first_min_quadobj.py $(fd pareto.txt results/migration_scan/$1) > results/migration_scan/$1/pareto/optima.tsv
