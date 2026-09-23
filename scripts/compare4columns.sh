#! /bin/bash

set -euf

trap 'rm -r $temp1' 0
trap 'rm -r $temp2' 0
temp1=$(mktemp)
temp2=$(mktemp)

cut -f 1-4 $1 > $temp1
cut -f 1-4 $2 > $temp2

diff $temp1 $temp2