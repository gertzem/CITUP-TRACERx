#! /bin/bash

set -euf

persons=( $(sed 1d "$1" | cut -f 1) )
dir=$2

for person in "${persons[@]}"; do
    echo "Reading $person"
    tree_scores=( $(find $dir -name "$person*.treeScores.tsv" | grep -v pareto) )
    if [ "${#tree_scores[@]}" -gt 0 ]; then
        mkdir -p "$dir/pareto/$person"
        python scripts/pareto_front.py "${tree_scores[@]}" > "$dir/pareto/$person/pareto.txt"
    fi
done
