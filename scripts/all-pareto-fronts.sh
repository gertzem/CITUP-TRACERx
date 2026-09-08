#! /bin/bash

set -euf

persons=( $(sed 1d "$1" | cut -f 1) )
dir=$2

for person in "${persons[@]}"; do
    echo "Reading $person"
    tree_scores=( $(find "$dir" -name "$person*.tree_0.txt" ! -name '*pareto*') )
    if [ "${#tree_scores[@]}" -gt 0 ]; then
        mkdir -p "$dir/pareto/$person"
        python scripts/read_citup_obj.py "${tree_scores[@]}" | sort -k 2,2n > "$dir/pareto/$person/pareto.txt"
    fi
done
