#! /bin/bash

set -uf
set -o pipefail

sourcedir=$1
( cd $1; python scripts/scan_termination_status.py tracerx ) | \
    grep -e TIMEOUT -e COMPLETED -e Infeasible | sort > a


while read person mig rest; do
    echo rm -rf results/migration_scan/tracerx/$mig/$person
    echo rm -rf results/migration_scan/tracerx/logs/$mig/$person
    echo rm -rf results/migration_scan/tracerx/pareto/$person
    echo rm -f results/pareto_plots/tracerx/${person}_pareto.pdf
    echo cp -a  "$1/results/migration_scan/tracerx/$mig/$person" "results/migration_scan/tracerx/$mig"
    echo cp -a  "$1/results/migration_scan/tracerx/logs/$mig/$person" "results/migration_scan/tracerx/logs/$mig"
done < a