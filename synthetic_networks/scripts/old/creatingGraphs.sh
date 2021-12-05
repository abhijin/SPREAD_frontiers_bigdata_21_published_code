#!/bin/bash
echo "scripting time"
prior=("python3" "maing.py" "--number_of_nodes" 144 --number_of_regions 9 --range 1 --locality_size 4 --locality_graph star --long_distance_type ER --ld_param .5)
p=("sbatch" "-o" "log" "--export=command=\"python3" "maing.py" "--number_of_nodes" 144 --number_of_regions 9 --range 1 --locality_size 4 --locality_graph star --long_distance_type ER --ld_param .5\" "run_proc.sbatch;")
localities=("star" "clique")
longDistanceTypes=("ER" "CL" "SF")
ERProbs=(0 0.25 0.5 0.75 1)
maxNumRegions=3 #rows of regions
numRegionsMax=$((maxNumRegions*maxNumRegions))
#for (( i=1; i <= $maxNumRegions; i++ ))
#do
#	p+=(0)
#done
#echo "${p[@]}"
for nodeRows in {2..2} #how many rows of nodes
do
	#for (( i=1; i <= $numRegionsMax; i++ ))
	#do
		#p[18+$i]=0
	#done
        for (( numberRegions=1; numberRegions <= $maxNumRegions; numberRegions++ )) #how many rows of regions
	do	
		#for (( i=1; i < $((numberRegions*numberRegions)); i++ ))
		#do
			#p[18+$i]=1 #need to adjust number when adding sbatch commands
		#done
		
                for range in {1..1} 
		do
                        for localityNodes in {1..2} #how many rows of locality nodes there are
			do
                                for localityType in {0..1}
				do
                                        for longDistance in {0..0}
					do
						for t in ${ERProbs[@]};
						do

	

						numNodes=$((nodeRows*nodeRows))
						localitySize=$((localityNodes*localityNodes))
						regionalNumber=$((numberRegions*numberRegions))
						#echo $numNodes
						if (( $numNodes>=regionalNumber )) && (( $(($numNodes/regionalNumber))>=$localitySize )) && (( $(($(($numNodes/regionalNumber))%2))==$(($localitySize%2)) ))
						then
							p[6]=$numNodes
							p[8]=$regionalNumber
							p[12]=$localitySize
							p[14]=${localities[$localityType]}
							p[16]=${longDistanceTypes[$longDistance]}
							p[18]=$t
							#"${p[@]}" #call sbatch command right here			
							#echo "${p[@]}"
sbatch -o log.${numNodes}.${regionalNumber}.${range}.${localitySize}.${localities[$localityType]}.${longDistanceTypes[$longDistance]}.$t --export=command="python3 ../scripts/maing.py --number_of_nodes $numNodes --number_of_regions $regionalNumber --range $range --locality_size $localitySize --locality_graph ${localities[$localityType]} --long_distance_type ${longDistanceTypes[$longDistance]} --ld_param $t" run_proc.sbatch;
qreg;
							#echo $numNodes
							#echo $regionalNumber
							echo ""
						fi
						done
                                        done
                                done	
                        done
                done
        done
done
#"${p[@]}" &> sum.txt
