#!/bin/bash
#SBATCH -c 24 # Number of cores reserved for each task
#SBATCH --mem=20G # Memory reserved for each task
#SBATCH --mail-type=END,FAIL # When to send an email
#SBATCH -J run_modwheel
#SBATCH --output=modwheel/run-%A.log
#SBATCH --constraint=titan-x
#SBATCH --time=72:00:00 # Time limit


HOSTNAME=`hostname`
NUM_CPUS=`lscpu | egrep '^CPU\(s\):' | cut -d ':' -f2- | sed -e 's/^[[:space:]]*//'`
NUM_CORES=`egrep '^cpu cores' /proc/cpuinfo | uniq | cut -d ':' -f2 | cut -d ' ' -f2-`
MEM=`lsmem | egrep 'Total online memory:' | cut -d: -f2- | sed -e 's/^[[:space:]]*//'`
echo "Run on ${HOSTNAME} which has ${NUM_CPUS} CPUs (${NUM_CORES} cores) and ${MEM} memory."

python3 run_modwheel.py 1 &
pid1=$!

python3 run_modwheel.py 2 &
pid2=$!

python3 run_modwheel.py 3 &
pid3=$!

python3 run_modwheel.py 4 &
pid4=$!

python3 run_modwheel.py 5 &
pid5=$!

python3 run_modwheel.py 6 &
pid6=$!

python3 run_modwheel.py 7 &
pid7=$!

python3 run_modwheel.py 8 &
pid8=$!

python3 run_modwheel.py 9 &
pid9=$!

python3 run_modwheel.py 10 &
pid10=$!

wait $pid1 $pid2 $pid3 $pid4 $pid5 $pid6 $pid7 $pid8 $pid9 $pid10

