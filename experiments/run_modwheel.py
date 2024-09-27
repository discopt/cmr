import sys
import os
import math
import subprocess
from datetime import datetime

from config_modwheel import *

assert os.path.exists(INSTANCE_DIRECTORY)

sample = int(sys.argv[1])
sampleStorage = f'{LOCAL_STORAGE}/modwheel_{sample}/'

os.system(f'mkdir -p {sampleStorage}')

def call(command):
#  print('[' + command + ']')
  os.system(command)

for instance in INSTANCES:
  file_base = f'{INSTANCE_DIRECTORY}/modwheel-{instance:05d}x{instance:05d}#{sample:03d}'

  print(f'Considering {file_base} at {datetime.now()}.', flush=True)

  # Run cmr-tu.
  call(f'gunzip -cd {file_base}.sparse.gz | {BUILD_DIRECTORY}/cmr-tu - -i sparse --stats --algo decomposition --time-limit 3600 1> {file_base}-cmrdec.out 2> {file_base}-cmrdec.err')
  if instance <= 50:
    call(f'gunzip -cd {file_base}.sparse.gz | {BUILD_DIRECTORY}/cmr-tu - -i sparse --stats --algo eulerian --time-limit 3600 1> {file_base}-cmreuler.out 2> {file_base}-cmreuler.err')
  if instance <= 50:
    call(f'gunzip -cd {file_base}.sparse.gz | {BUILD_DIRECTORY}/cmr-tu - -i sparse --stats --algo partition --time-limit 3600 1> {file_base}-cmrpart.out 2> {file_base}-cmrpart.err')
    
  call(f'gunzip -cd {file_base}.sparse.gz | {BUILD_DIRECTORY}/cmr-tu - -i sparse --stats --algo decomposition --time-limit 3600 -N {file_base}-cmrcert.sub 1> {file_base}-cmrcert.out 2> {file_base}-cmrcert.err')

  # Run unimodularity-test.
  call(f'gunzip -cd {file_base}.sparse.gz | {BUILD_DIRECTORY}/cmr-matrix - -i sparse -o dense {sampleStorage}/input.dense')
  call(f'{UNIMOD_DIRECTORY}/unimodularity-test {sampleStorage}/input.dense -s 2> /dev/null | egrep \'^[ 0-9-]*$\' 1> {sampleStorage}/signed.dense')
  call(f'{UNIMOD_DIRECTORY}/unimodularity-test {sampleStorage}/signed.dense -t -v 1> {file_base}-unimod.out 2> {file_base}-unimod.err')
  if instance <= 199:
    call(f'{UNIMOD_DIRECTORY}/unimodularity-test {sampleStorage}/signed.dense -t -v -c 1> {file_base}-unimodcert.out 2> {file_base}-unimodcert.err')

os.system(f'rm -r {sampleStorage}')

