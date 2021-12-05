#!/bin/bash
# SLURM JOB SCRIPT OPTIONS:
#SBATCH --time=24:00:00

##SBATCH --mem=30464
##SBATCH --mem-per-cpu=1904


# Set the partition to submit to (a partition is equivalent to a queue)
#SBATCH -p bii
#SBATCH --qos bii-200

# charge code
#SBATCH -A bii_nssac
##SBATCH -A nssac_students

# Load any modules you require:

echo "###################################"
echo "$command"
echo "###################################"
eval "$command"
