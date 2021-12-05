#!/bin/bash
function test(){
python ../scripts/generate_instances.py ../input/test_design.json
}

function real(){
## python ../scripts/generate_instances.py ../input/BD.json
## python ../scripts/generate_instances.py ../input/ID.json
## python ../scripts/generate_instances.py ../input/TH.json
python ../scripts/generate_instances.py ../input/VN.json
}

function design(){
for s in `seq 1 10`
do
    for l in `echo 4`
    do
        for p in `echo 0.1`
        do
            for seed in `echo all`
            do
            cat << EOF > s_${s}_seed_${seed}_design.json
{
    "network": "/sfs/qumulo/qproject/biocomplexity/abhijin/SPREAD_frontiers_bigdata/synthetic_networks/work/k_16/r_4/s_${l}/HL_clique/FLD_ER/p_${p}/seed_${s}",
    "unwanted_network_prefix": "/sfs/qumulo/qproject/biocomplexity/abhijin/SPREAD_frontiers_bigdata/synthetic_networks/work/",
    "time_steps": 24,
    "dag_type": 0,
    "seeding_specification": "folder",

    "experiments": {
        "exposure_delay": [0], 
        "alpha_S": [0, 0.001, 0.005, 0.01, 0.015],
        "alpha_L": [0, 0.001, 0.005, 0.01, 0.015],
        "alpha_LD": [0, 0.001, 0.005, 0.01, 0.015],
        "number_of_simulations": [100],
        "moore_range": [1, 2, 3],
        "start_month": [1],
        "random_seed": [1],
        "seeding": "/sfs/qumulo/qproject/biocomplexity/abhijin/SPREAD_frontiers_bigdata/synthetic_networks/work/k_16/r_4/s_${l}/HL_clique/FLD_ER/p_${p}/seed_${s}/seed_${seed}.csv"
    }
}
EOF
            done
        done
    done
done
}

function run(){
for f in `find ../input/net_16_4_p0.1/ -iname 's_*json'`
do
fname=`echo $f | sed -e 's/.*_p/p/' -e 's|/|-|g'`
echo $fname
python ../scripts/generate_instances.py $f
mv run.sh run.sh_$fname
done
cat run.sh_* > run.sh
rm run.sh_*
}

function collect_results() {    # collect results from log files
find $EXPERIMENTS -iname "log" | xargs -I {} grep '^INSERT' {} > to_db.sqlite
find $EXPERIMENTS -iname "log" | xargs grep -L '^INSERT' | xargs -I {} tail -n1 -v {} > to_be_run_again.txt
echo "$(cat << EOF 
Check to_db.sqlite. If everything seems okay, then do the following ...
sqlite3 $RESULTS_DB < to_db.sqlite
Runs which did not go through are in to_be_run_again.txt. Run 'master.sh generate_instances' to rerun them.
EOF
)"
wc -l to_db.sqlite
wc -l to_be_run_again.txt
}

function import_to_db() {   # import results to DB
sqlite3 $RESULTS_DB < to_db.sqlite
}

function create_tables() {  # One time only: create results database
sqlite3 $RESULTS_DB << 'END_SQL'
CREATE TABLE IF NOT EXISTS "summary" (
    network TEXT,
    random_seed INTEGER,
    suitability_thresh REAL,
    exposure_delay INTEGER,
    moore_range INTEGER,
    alpha_S REAL,
    alpha_L REAL,
    alpha_LD REAL,
    start_month INTEGER,
    number_of_time_steps INTEGER,
    number_of_simulations INTEGER,
    seeding TEXT,
    interventions TEXT,
    interventions_type TEXT,
    time_step INTEGER,
    accumulated_probabilities REAL,
    infections_mean REAL,
    infections_std REAL,
    infections_min REAL,
    infections_25_per REAL,
    infections_50_per REAL,
    infections_75_per REAL,
    infections_max REAL);
END_SQL
##    PRIMARY KEY(network,random_seed,suitability_thresh,exposure_delay,moore_range,alpha_S,alpha_L,alpha_LD,start_month,number_of_time_steps,number_of_simulations,seeding,interventions,interventions_type,time_step));
}

function create_tables_exhaustive() {  # One time only: create results database exhaustive results
sqlite3 $RESULTS_DB_E << 'END_SQL'
CREATE TABLE IF NOT EXISTS "summary" (
    network TEXT,
    random_seed INTEGER,
    suitability_thresh REAL,
    exposure_delay INTEGER,
    moore_range INTEGER,
    alpha_S REAL,
    alpha_L REAL,
    alpha_LD REAL,
    start_month INTEGER,
    number_of_time_steps INTEGER,
    number_of_simulations INTEGER,
    seeding TEXT,
    interventions TEXT,
    interventions_type TEXT,
    time_step INTEGER,
    accumulated_probabilities REAL,
    infections_mean REAL,
    infections_std REAL,
    infections_min REAL,
    infections_25_per REAL,
    infections_50_per REAL,
    infections_75_per REAL,
    infections_max REAL);
END_SQL
}

function import_to_db_E() {   # import results to DB
sqlite3 $RESULTS_DB_E < to_db.sqlite
}

if [[ $# == 0 ]]; then
    echo "Assumes the following environment variables are set:"
    echo "SIMULATOR; EXPERIMENTS; RESULTS_DB"
    echo "Here are the options:"
    grep "^function" $BASH_SOURCE | sed -e 's/function/  /' -e 's/[(){]//g' -e '/IGNORE/d'
else
    if [ "$SIMULATOR" = "" ] || [ "$EXPERIMENTS" = "" ] || [ "$RESULTS_DB" = "" ] ; then
        echo "ERROR: Enviroment variables not set. See scripts/setenv/."
        exit
    fi
    eval $1 $2
fi
