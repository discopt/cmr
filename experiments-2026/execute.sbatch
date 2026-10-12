#!/bin/bash
#SBATCH -c 24 # Number of cores reserved for each task
#SBATCH --mem=120G # Memory reserved for each task
#SBATCH --mail-type=END,FAIL # When to send an email
#SBATCH -J cmr
#SBATCH --output=results/main-%A.log
#SBATCH --constraint=titan-x
#SBATCH --time=168:00:00 # Time limit


HOSTNAME=`hostname`
NUM_CPUS=`lscpu | egrep '^CPU\(s\):' | cut -d ':' -f2- | sed -e 's/^[[:space:]]*//'`
NUM_CORES=`egrep '^cpu cores' /proc/cpuinfo | uniq | cut -d ':' -f2 | cut -d ' ' -f2-`
MEM=`lsmem | egrep 'Total online memory:' | cut -d: -f2- | sed -e 's/^[[:space:]]*//'`
UPTIME=`uptime`
echo "Run on ${HOSTNAME} which has ${NUM_CPUS} CPUs (${NUM_CORES} cores) and ${MEM} memory. Uptime: ${UPTIME}"

# 4GB
ulimit -Sv 12000000

python3 execute.py &
pid1=$!

python3 execute.py &
pid2=$!

python3 execute.py &
pid3=$!

python3 execute.py &
pid4=$!

python3 execute.py &
pid5=$!

python3 execute.py &
pid6=$!

python3 execute.py &
pid7=$!

python3 execute.py &
pid8=$!

python3 execute.py &
pid9=$!

python3 execute.py &
pid10=$!

wait $pid1 $pid2 $pid3 $pid4 $pid5 $pid6 $pid7 $pid8 $pid9 $pid10

UPTIME=`uptime`
echo "Ran on ${HOSTNAME} which has ${NUM_CPUS} CPUs (${NUM_CORES} cores) and ${MEM} memory. Uptime: ${UPTIME}"
