#!/bin/bash
#SBATCH -c 1 # Number of cores reserved for each task
#SBATCH --mem=1G # Memory reserved for each task
#SBATCH --mail-type=END,FAIL # When to send an email
#SBATCH -J gen_network
#SBATCH --output=network/gen-%A.log
#SBATCH --time=8:00:00 # Time limit


HOSTNAME=`hostname`
NUM_CPUS=`lscpu | egrep '^CPU\(s\):' | cut -d ':' -f2- | sed -e 's/^[[:space:]]*//'`
NUM_CORES=`egrep '^cpu cores' /proc/cpuinfo | uniq | cut -d ':' -f2 | cut -d ' ' -f2-`
MEM=`lsmem | egrep 'Total online memory:' | cut -d: -f2- | sed -e 's/^[[:space:]]*//'`
echo "Run on ${HOSTNAME} which has ${NUM_CPUS} CPUs (${NUM_CORES} cores) and ${MEM} memory."

python3 gen_network.py 10
