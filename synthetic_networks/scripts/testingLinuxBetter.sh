#!/usr/bin/
echo "scripting time"
p=("python3" "maing.py" "--number_of_nodes" 144 --number_of_regions 9 --range 1 --locality_size 4 --locality_graph star --long_distance_type ER --ld_param .5)
localities=("star" "clique")
longDistanceTypes=("ER" "CL" "SF")
ERProbs=(0 0.25 0.5 0.75 1)
maxNumRegions=3 #rows of regions
numRegionsMax=$((maxNumRegions*maxNumRegions))
for (( i=1; i <= $maxNumRegions; i++ ))
do
	p+=(0)
done
echo "${p[@]}"
for nodeRows in {2..3} #how many rows of nodes
do
	for (( i=1; i <= $numRegionsMax; i++ ))
	do
		p[15+$i]=0
	done
        for (( numberRegions=1; numberRegions <= $maxNumRegions; numberRegions++ )) #how many rows of regions
	do	
		for (( i=1; i < $((numberRegions*numberRegions)); i++ ))
		do
			p[15+$i]=1 #need to adjust number when adding sbatch commands
		done
		
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
							p[3]=$numNodes
							p[5]=$regionalNumber
							p[9]=$localitySize
							p[11]=${localities[$localityType]}
							p[13]=${longDistanceTypes[$longDistance]}
							p[15]=$t
							"${p[@]}" #call sbatch command right here			
							#echo "${p[@]}"
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