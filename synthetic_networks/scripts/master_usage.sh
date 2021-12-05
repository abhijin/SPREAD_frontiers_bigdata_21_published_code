#!/bin/bash
function netgen(){
python3 ../scripts/mpsn.py \
--number_of_nodes 4096 \
--number_of_regions 16 \
--range 3 \
--locality_size 4 \
--locality_graph clique \
--long_distance_type ER \
--seed 0 \
--ld_param 1
}

if [[ $# == 0 ]]; then
    grep "^function" $BASH_SOURCE | sed -e 's/function/  /' -e 's/[(){]//g' -e '/IGNORE/d'
    exit
fi
eval $1 $2
