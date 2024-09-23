import sys
import os
import math
import subprocess

from config import *
INSTANCE_DIRECTORY = 'network'

assert os.path.exists(INSTANCE_DIRECTORY)

os.system(f'mkdir -p {LOCAL_STORAGE}/network')

def call(command):
#  print('[' + command + ']')
  os.system(command)

for order in sorted(list(range(1,41)) + [ 100 * i for i in range(1,41) ] + [ 1000 * i for i in range(1,41) ]):
  r = 1
  while True:
    file_base = f'{INSTANCE_DIRECTORY}/network-{order:05d}x{order:05d}#{r:03d}'
    if not os.path.exists(f'{file_base}.sparse.gz'):
        break

    print(f'Considering {file_base}.')

    # Run cmr-tu.
    call(f'gunzip -cd {file_base}.sparse.gz | {BUILD_DIRECTORY}/cmr-tu - -i sparse --stats --algo decomposition --time-limit 3600 1> {file_base}-cmrdec.out 2> {file_base}-cmrdec.err')
    if order <= 20:
      call(f'gunzip -cd {file_base}.sparse.gz | {BUILD_DIRECTORY}/cmr-tu - -i sparse --stats --algo eulerian --time-limit 3600 1> {file_base}-cmreuler.out 2> {file_base}-cmreuler.err')
    if order <= 30:
      call(f'gunzip -cd {file_base}.sparse.gz | {BUILD_DIRECTORY}/cmr-tu - -i sparse --stats --algo partition --time-limit 3600 1> {file_base}-cmrpart.out 2> {file_base}-cmrpart.err')

    # Run unimodularity-test.
    if order <= 4000:
      call(f'gunzip -cd {file_base}.sparse.gz | {BUILD_DIRECTORY}/cmr-matrix - -i sparse -o dense {LOCAL_STORAGE}/network/input.dense')
      call(f'{UNIMOD_DIRECTORY}/unimodularity-test {LOCAL_STORAGE}/network/input.dense -t -v 1> {file_base}-unimod.out 2> {file_base}-unimod.err')

    r += 1

os.system(f'rm -r {LOCAL_STORAGE}/network')

