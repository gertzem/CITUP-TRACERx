#! /bin/bash

set -euf

logs=( $(find $1 -name 'citup.log' | grep -v pareto) )
for ell in ${logs[@]}; do
    d=$(dirname $ell)
    p=$(basename $d)
    rt=$(dirname $(dirname $(dirname $d)))
    m=$(basename $(dirname $d))
    if ! compgen -G $rt/$m/$p/*.treeScores.tsv > /dev/null; then
        if grep -q 'odel is infeasible' $ell; then
            echo $p $m infeasible
        else
            echo $p $m running
        fi
    fi
done