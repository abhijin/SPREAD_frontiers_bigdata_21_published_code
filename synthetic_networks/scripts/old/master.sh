#!/bin/bash
function mpn(){ # usage of multipathway network script
    python ../scripts/multipathway_net.py \
        --number_of_nodes 64 \
        --number_of_regions 4 \
        --range 1 \
        --locality_size 4 \
        --locality_graph clique \
        --long_distance_type ER \
        --ld_param .5
}


if [[ $# == 0 ]]; then
   echo "Here are the options:"
   grep "^function" $BASH_SOURCE | sed -e 's/function/  /' -e 's/[(){]//g' -e '/IGNORE/d'
else
   eval $1 $2
fi

