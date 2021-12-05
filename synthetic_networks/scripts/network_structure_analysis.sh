#!/bin/bash
number_of_rows=64
range_list="1 1.5 2 2.5 3 3.5 4"
regions_per_row_list="2 4 8 16"
locality_rows_list="2 4 8"
locality_graph_list="star clique"
long_distance_graph_list="ER" # CL SF"
long_distance_epsilon_list="0.1 .5 1 5 10"
script=`realpath ../scripts/mpsn.py`

((count=0))
for range in `echo "$range_list"`
do
    for regions_per_row in `echo "$regions_per_row_list"`
    do	
        for locality_rows in `echo "$locality_rows_list"`
    	do
            for locality_graph in `echo "$locality_graph_list"`
    		do
                for long_distance_graph in `echo "$long_distance_graph_list"`
    			do
        			for eps in `echo "$long_distance_epsilon_list"`
    				do
                        for seed in `seq 1 10`
                        do
                            ((count++))
                            number_of_nodes=$((number_of_rows*number_of_rows))
    				        number_of_regions=$((regions_per_row*regions_per_row))
    				        locality_size=$((locality_rows*locality_rows))
                            folder=k_${number_of_regions}/r_${range}/s_${locality_size}/HL_${locality_graph}/FLD_${long_distance_graph}/p_$eps/seed_$seed

                            skip_flag=`grep -l INSERT $folder/log | wc -l`
                            if [[ $skip_flag == 1 ]]; then
                                continue
                            fi
                            mkdir -p $folder
disp=`sbatch -o log -D $folder --export=command="python3 $script \
--number_of_nodes $number_of_nodes \
--number_of_regions $number_of_regions \
--range $range \
--locality_size $locality_size \
--locality_graph $locality_graph \
--long_distance_type $long_distance_graph \
--seed $seed \
--ld_param $eps" ../scripts/run_proc.sbatch;`
                        echo $count $disp
bash ../scripts/qreg;
                        done
                    done
                done
            done
        done
    done
done
