#!/bin/bash
function eig_diam(){
    python ../scripts/eig_diam.py
}

function unravel(){
python ../scripts/unravel_pathways.py --range 2 --locality_size 16 --epsilon 1

python ../scripts/unravel_pathways.py --range 1 --locality_size 16 --epsilon 1
python ../scripts/unravel_pathways.py --range 3 --locality_size 16 --epsilon 1

python ../scripts/unravel_pathways.py --range 2 --locality_size 16 --epsilon 0.1
python ../scripts/unravel_pathways.py --range 2 --locality_size 16 --epsilon 10

python ../scripts/unravel_pathways.py --range 2 --locality_size 4 --epsilon 1
}

function non_uniform(){
python ../scripts/non_uniform_pathways.py --range 2 --locality_size 16 --epsilon 1
python ../scripts/non_uniform_pathways.py --range 2 --locality_size 16 --epsilon 0.1
python ../scripts/non_uniform_pathways.py --range 2 --locality_size 16 --epsilon 10
python ../scripts/non_uniform_pathways.py --range 1 --locality_size 16 --epsilon 1
python ../scripts/non_uniform_pathways.py --range 3 --locality_size 16 --epsilon 1
python ../scripts/non_uniform_pathways.py --range 2 --locality_size 4 --epsilon 1
}

function transfer(){
    mv unravel_s_16_eps1.0_range1.pdf ../../../frontiers_bigdata_ias/figs/
    mv unravel_s_16_eps1.0_range2.pdf ../../../frontiers_bigdata_ias/figs/
    mv unravel_s_16_eps1.0_range3.pdf ../../../frontiers_bigdata_ias/figs/
    mv unravel_s_16_eps10.0_range2.pdf ../../../frontiers_bigdata_ias/figs/
    mv unravel_s_16_eps0.1_range2.pdf ../../../frontiers_bigdata_ias/figs/
    mv unravel_s_4_eps1.0_range2.pdf ../../../frontiers_bigdata_ias/figs/
    mv non_uniform_s_16_eps0.1_range2.pdf  ../../../frontiers_bigdata_ias/figs/
    mv non_uniform_s_16_eps1.0_range1.pdf  ../../../frontiers_bigdata_ias/figs/
    mv non_uniform_s_16_eps1.0_range2.pdf  ../../../frontiers_bigdata_ias/figs/
    mv non_uniform_s_16_eps1.0_range3.pdf  ../../../frontiers_bigdata_ias/figs/
    mv non_uniform_s_16_eps10.0_range2.pdf  ../../../frontiers_bigdata_ias/figs/
    mv non_uniform_s_4_eps1.0_range2.pdf  ../../../frontiers_bigdata_ias/figs/
    mv real_diameter_1.pdf ../../../frontiers_bigdata_ias/figs/
    mv real_diameter_2.pdf ../../../frontiers_bigdata_ias/figs/
    mv real_unweighted_spectral_radius_1.pdf ../../../frontiers_bigdata_ias/figs/
    mv real_unweighted_spectral_radius_2.pdf ../../../frontiers_bigdata_ias/figs/
    mv real_weighted_spectral_radius_1.pdf ../../../frontiers_bigdata_ias/figs/
    mv real_weighted_spectral_radius_2.pdf ../../../frontiers_bigdata_ias/figs/
}

if [[ $# == 0 ]]; then
    echo "Here are the options:"
    grep "^function" $BASH_SOURCE | sed -e 's/function/  /' -e 's/[(){]//g' -e '/IGNORE/d'
else
    eval $1 $2
fi
