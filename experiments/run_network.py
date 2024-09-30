import sys
import os
import math
import subprocess

from config_network import *

assert os.path.exists(INSTANCE_DIRECTORY)

sample = int(sys.argv[1])
sampleStorage = f'{LOCAL_STORAGE}/network_{sample}/'

def call(command):
#  print('[' + command + ']')
  os.system(command)

for instance in INSTANCES:
  file_base = f'{INSTANCE_DIRECTORY}/network-{order:05d}x{order:05d}#{sample:03d}'

  print(f'Considering {file_base}.')

  # Run cmr-tu.
  call(f'gunzip -cd {file_base}.sparse.gz | {BUILD_DIRECTORY}/cmr-tu - -i sparse --stats --algo decomposition --time-limit 3600 1> {file_base}-cmrdec.out 2> {file_base}-cmrdec.err')
  if order <= 20:
    call(f'gunzip -cd {file_base}.sparse.gz | {BUILD_DIRECTORY}/cmr-tu - -i sparse --stats --algo eulerian --time-limit 3600 1> {file_base}-cmreuler.out 2> {file_base}-cmreuler.err')
  if order <= 30:
    call(f'gunzip -cd {file_base}.sparse.gz | {BUILD_DIRECTORY}/cmr-tu - -i sparse --stats --algo partition --time-limit 3600 1> {file_base}-cmrpart.out 2> {file_base}-cmrpart.err')

  # Run unimodularity-test.
  if order <= 4000:
    call(f'gunzip -cd {file_base}.sparse.gz | {BUILD_DIRECTORY}/cmr-matrix - -i sparse -o dense {sampleStorage}/input.dense')
    call(f'{UNIMOD_DIRECTORY}/unimodularity-test {sampleStorage}/input.dense -t -v 1> {file_base}-unimod.out 2> {file_base}-unimod.err')

os.system(f'rm -r {sampleStorage}')

