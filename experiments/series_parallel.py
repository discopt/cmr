import sys
import os
import math
import subprocess

GENERATOR = '../build-release/cmr-generate-series-parallel'

# Parameter of script: max amount of GB to be used.
try:
  maxMemory = float(sys.argv[1]) * 1024 * 1024 * 1024
  numRepetitions = int(sys.argv[2])
except:
  print(f'Usage: {sys.argv[0]} MAX-MEMORY-IN-GIGABYTES NUM-REPETITIONS')
  sys.exit(1)
maxNumNonzeros = maxMemory / 45 # This is the approximate number of bytes required per nonzero.

bitsNonzeros = int(round(math.log(maxNumNonzeros) / math.log(2), 0))

sys.stdout.write('rTot,cTot,rBase,cBase,rZero,cZero,rUnit,cUnit,rCopy,cCopy,sp,tRed,tWheel,tTern,tTotal,#nzs\n')
sys.stdout.flush()

def run(sparsity, numBaseRows, numBaseColumns, numZeroRows, numZeroColumns, numUnitRows, numUnitColumns, numCopiedRows, numCopiedColumns, ternary, repetitions):
  command = [GENERATOR, str(int(numBaseRows)), str(int(numBaseColumns)), '-z', str(int(numZeroRows)), str(int(numZeroColumns)), '-u', str(int(numUnitRows)), str(int(numUnitColumns)), '-c', str(int(numCopiedRows)), str(int(numCopiedColumns)), '-s', str(sparsity), '-b', str(repetitions)]
  if ternary:
    command.append('-t')
  process = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
  timeReductions = 0.0
  timeWheel = 0.0
  timeTernary = 0.0
  timeTotal = 0.0
  avgNonzeros = None
  for line in process.stdout.decode('utf-8').split('\n'):
    if line.startswith('Search for reductions:'):
      timeReductions = float(line.split('/')[1].strip())
    if line.startswith('Search for wheel matrices:'):
      timeWheel = float(line.split('/')[1].strip())
    if line.startswith('Search for ternary certificate:'):
      timeTernary = float(line.split('/')[1].strip())
    if line.startswith('Total:'):
      timeTotal = float(line.split('/')[1].strip())
    if line.startswith('Average number of nonzeros:'):
      avgNonzeros = float(line.split(':')[1].strip())
  timeReductions /= repetitions
  timeWheel /= repetitions
  timeTernary /= repetitions
  timeTotal /= repetitions
  sys.stdout.write(f'{numBaseRows+numZeroRows+numUnitRows+numCopiedRows:.0f},{numBaseColumns+numZeroColumns+numUnitColumns+numCopiedColumns:.0f},{numBaseRows},{numBaseColumns},{numZeroRows},{numZeroColumns},{numUnitRows},{numUnitColumns},{numCopiedRows},{numCopiedColumns},{sparsity},{timeReductions},{timeWheel},{timeTernary},{timeTotal},{avgNonzeros}\n')
  sys.stdout.flush()

for ternary in [False, True]:

  # Different portions of unit/copied.
  if True:
    size = 2**(bitsNonzeros-16)
    run(1, 1, 1, 0, 0, 1.0*size-1, 1.0*size-1, 0, 0, ternary, numRepetitions)
    run(1, 1, 1, 0, 0, 0.9*size-1, 0.9*size-1, 0.1*size, 0.1*size, ternary, numRepetitions)
    run(1, 1, 1, 0, 0, 0.8*size-1, 0.8*size-1, 0.2*size, 0.2*size, ternary, numRepetitions)
#run(1, 1, 1, 0, 0, 0.7*size-1, 0.7*size-1, 0.3*size, 0.3*size, ternary, numRepetitions)
#run(1, 1, 1, 0, 0, 0.6*size-1, 0.6*size-1, 0.4*size, 0.4*size, ternary, numRepetitions)
#run(1, 1, 1, 0, 0, 0.5*size-1, 0.5*size-1, 0.5*size, 0.5*size, ternary, numRepetitions)
#run(1, 1, 1, 0, 0, 0.4*size-1, 0.4*size-1, 0.6*size, 0.6*size, ternary, numRepetitions)
#run(1, 1, 1, 0, 0, 0.3*size-1, 0.3*size-1, 0.7*size, 0.7*size, ternary, numRepetitions)
    run(1, 1, 1, 0, 0, 0.2*size-1, 0.2*size-1, 0.8*size, 0.8*size, ternary, numRepetitions)
    run(1, 1, 1, 0, 0, 0.1*size-1, 0.1*size-1, 0.9*size, 0.9*size, ternary, numRepetitions)
    run(1, 1, 1, 0, 0, 0, 0, 1.0*size-1, 1.0*size-1, ternary, numRepetitions)
    sys.stdout.write('\n')
  
  if False:
    # Different layouts, but same number of nonzeros.
    size = 2**(bitsNonzeros-16)
    run(1, 1, 1, 0, 0, 0.5*size/8, 0.5*size*2.3, 0.5*size/8, 0.5*size*2.3, ternary, numRepetitions)
    run(1, 1, 1, 0, 0, 0.5*size/4, 0.5*size*2.0, 0.5*size/4, 0.5*size*2.0, ternary, numRepetitions)
    run(1, 1, 1, 0, 0, 0.5*size/2, 0.5*size*1.6, 0.5*size/2, 0.5*size*1.6, ternary, numRepetitions)
    run(1, 1, 1, 0, 0, 0.5*size, 0.5*size, 0.5*size, 0.5*size, ternary, numRepetitions)
    run(1, 1, 1, 0, 0, 0.5*size*1.6, 0.5*size/2, 0.5*size*1.6, 0.5*size/2, ternary, numRepetitions)
    run(1, 1, 1, 0, 0, 0.5*size*2.0, 0.5*size/4, 0.5*size*2.0, 0.5*size/4, ternary, numRepetitions)
    run(1, 1, 1, 0, 0, 0.5*size*2.3, 0.5*size/8, 0.5*size*2.3, 0.5*size/8, ternary, numRepetitions)
    sys.stdout.write('\n')
  
  if False:
    # Different portion of base.
    size = 2**(bitsNonzeros-16)
    run(1, 1, 1, 0, 0, 0.5*size-1, 0.5*size-1, 0.5*size, 0.5*size, ternary, numRepetitions)
    run(0.05*size, 0.1*size, 0.1*size, 0, 0, 0.45*size, 0.45*size, 0.45*size, 0.45*size, ternary, numRepetitions)
    run(0.10*size, 0.2*size, 0.2*size, 0, 0, 0.40*size, 0.40*size, 0.40*size, 0.40*size, ternary, numRepetitions)
    run(0.15*size, 0.3*size, 0.3*size, 0, 0, 0.35*size, 0.35*size, 0.35*size, 0.35*size, ternary, numRepetitions)
    run(0.20*size, 0.4*size, 0.4*size, 0, 0, 0.30*size, 0.30*size, 0.30*size, 0.30*size, ternary, numRepetitions)
    run(0.25*size, 0.5*size, 0.5*size, 0, 0, 0.25*size, 0.25*size, 0.25*size, 0.25*size, ternary, numRepetitions)
    run(0.30*size, 0.6*size, 0.6*size, 0, 0, 0.20*size, 0.20*size, 0.20*size, 0.20*size, ternary, numRepetitions)
    run(0.35*size, 0.7*size, 0.7*size, 0, 0, 0.15*size, 0.15*size, 0.15*size, 0.15*size, ternary, numRepetitions)
    run(0.40*size, 0.8*size, 0.8*size, 0, 0, 0.10*size, 0.10*size, 0.10*size, 0.10*size, ternary, numRepetitions)
    run(0.45*size, 0.9*size, 0.9*size, 0, 0, 0.05*size, 0.05*size, 0.05*size, 0.05*size, ternary, numRepetitions)
    run(0.50*size, 1.0*size, 1.0*size, 0, 0, 0, 0, 0, 0, ternary, numRepetitions)
    sys.stdout.write('\n')
    
  if False:
    # Different densities of full base.
    size = 2**(bitsNonzeros-16)
    run(0.0*size, size, size, 0, 0, 0, 0, 0, 0, ternary, numRepetitions)
    run(0.1*size, size, size, 0, 0, 0, 0, 0, 0, ternary, numRepetitions)
    run(0.2*size, size, size, 0, 0, 0, 0, 0, 0, ternary, numRepetitions)
    run(0.3*size, size, size, 0, 0, 0, 0, 0, 0, ternary, numRepetitions)
    run(0.4*size, size, size, 0, 0, 0, 0, 0, 0, ternary, numRepetitions)
    run(0.5*size, size, size, 0, 0, 0, 0, 0, 0, ternary, numRepetitions)
    run(0.6*size, size, size, 0, 0, 0, 0, 0, 0, ternary, numRepetitions)
    run(0.7*size, size, size, 0, 0, 0, 0, 0, 0, ternary, numRepetitions)
    run(0.8*size, size, size, 0, 0, 0, 0, 0, 0, ternary, numRepetitions)
    run(0.9*size, size, size, 0, 0, 0, 0, 0, 0, ternary, numRepetitions)
    run(1.0*size, size, size, 0, 0, 0, 0, 0, 0, ternary, numRepetitions)
    sys.stdout.write('\n')
    
