#!/bin/bash
#SBATCH -c 1 # Number of cores reserved for each task
#SBATCH --mem=40G # Memory reserved for each task
#SBATCH --mail-type=END,FAIL # When to send an email
#SBATCH -J run_network
#SBATCH --output=network/run-%A.log
#SBATCH --constraint=titan-x
#SBATCH --time=168:00:00 # Time limit


HOSTNAME=`hostname`
NUM_CPUS=`lscpu | egrep '^CPU\(s\):' | cut -d ':' -f2- | sed -e 's/^[[:space:]]*//'`
NUM_CORES=`egrep '^cpu cores' /proc/cpuinfo | uniq | cut -d ':' -f2 | cut -d ' ' -f2-`
MEM=`lsmem | egrep 'Total online memory:' | cut -d: -f2- | sed -e 's/^[[:space:]]*//'`
echo "Run on ${HOSTNAME} which has ${NUM_CPUS} CPUs (${NUM_CORES} cores) and ${MEM} memory."

for SAMPLE in `seq 1 10`; do
  python3 run_network.py ${SAMPLE}
done
