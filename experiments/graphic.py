import sys
import os
import math
import subprocess

GENERATOR = '../build-release/cmr-generate-graphic'

# Parameter of script: #rows and number of repetitions
try:
  numRows = int(sys.argv[1])
  numRepetitions = int(sys.argv[2])
except:
  print(f'Usage: {sys.argv[0]} #ROWS NUM-REPETITIONS')
  sys.exit(1)

sys.stdout.write('rows,cols,#nzs,tTrans,tCheck,tApply,tTotal\n')
sys.stdout.flush()

def run(numRows, numColumns, repetitions):
  command = [GENERATOR, str(int(numRows)), str(int(numColumns)), '-b', str(repetitions)]
  process = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
  timeTranspose = 0.0
  timeCheck = 0.0
  timeApply = 0.0
  timeTotal = 0.0
  avgNonzeros = 0
  for line in process.stdout.decode('utf-8').split('\n'):
    if line.startswith('Generated a'):
      avgNonzeros += int(line.split(' ')[5].strip())
    if line.startswith('Transposition:'):
      timeTranspose += float(line.split('/')[1].strip())
    if line.startswith('Check:'):
      timeCheck += float(line.split('/')[1].strip())
    if line.startswith('Apply:'):
      timeApply += float(line.split('/')[1].strip())
    if line.startswith('Total:'):
      timeTotal += float(line.split('/')[1].strip())
  timeTranspose /= repetitions
  timeCheck /= repetitions
  timeApply /= repetitions
  timeTotal /= repetitions
  avgNonzeros /= repetitions
  sys.stdout.write(f'{numRows:.0f},{numColumns:.0f},{avgNonzeros:.1f},{timeTranspose},{timeCheck},{timeApply},{timeTotal}\n')
  sys.stdout.flush()

run(numRows, int(0.5*numRows), numRepetitions)
run(numRows, int(1.0*numRows), numRepetitions)
run(numRows, int(1.5*numRows), numRepetitions)
run(numRows, int(2.0*numRows), numRepetitions)
run(numRows, int(2.5*numRows), numRepetitions)
run(numRows, int(3.0*numRows), numRepetitions)
run(numRows, int(3.5*numRows), numRepetitions)
run(numRows, int(4.5*numRows), numRepetitions)
run(numRows, int(5.0*numRows), numRepetitions)
run(numRows, int(5.5*numRows), numRepetitions)
run(numRows, int(6.0*numRows), numRepetitions)
run(numRows, int(6.5*numRows), numRepetitions)
run(numRows, int(7.0*numRows), numRepetitions)
run(numRows, int(7.5*numRows), numRepetitions)
run(numRows, int(8.0*numRows), numRepetitions)

