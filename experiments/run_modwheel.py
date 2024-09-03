import sys
import os
import math
import subprocess

INSTANCE_DIRECTORY = 'modwheel'
BUILD_DIRECTORY = '../build-release'
UNIMOD_DIRECTORY = '/home/matthias/work/cmr/unimodularity-test/unimodularity-library-1.2h/src/'

assert os.path.exists(INSTANCE_DIRECTORY)

def call(command):
#  print('[' + command + ']')
  os.system(command)

for order in sorted(list(range(3, 50, 2)) + list(range(99, 1000, 50))):
  r = 1
  while True:
    file_base = f'{INSTANCE_DIRECTORY}/modwheel-{order:05d}x{order:05d}#{r:03d}'
    print(f'{file_base}.sparse.gz')
    if not os.path.exists(f'{file_base}.sparse.gz'):
        break


    # Run cmr-tu.
    call(f'gunzip -cd {file_base}.sparse.gz | {BUILD_DIRECTORY}/cmr-tu - -i sparse --stats --algo decomposition --time-limit 3600 1> {file_base}-cmrdec.out 2> {file_base}-cmrdec.err')
    if order <= 40:
      call(f'gunzip -cd {file_base}.sparse.gz | {BUILD_DIRECTORY}/cmr-tu - -i sparse --stats --algo eulerian --time-limit 3600 1> {file_base}-cmreuler.out 2> {file_base}-cmreuler.err')
    if order <= 50:
      call(f'gunzip -cd {file_base}.sparse.gz | {BUILD_DIRECTORY}/cmr-tu - -i sparse --stats --algo partition --time-limit 3600 1> {file_base}-cmrpart.out 2> {file_base}-cmrpart.err')
    
    call(f'gunzip -cd {file_base}.sparse.gz | {BUILD_DIRECTORY}/cmr-tu - -i sparse --stats --algo decomposition --time-limit 3600 -N {file_base}-cmrcert.sub 1> {file_base}-cmrcert.out 2> {file_base}-cmrcert.err')

    # Run unimodularity-test.
    call(f'gunzip -cd {file_base}.sparse.gz | {BUILD_DIRECTORY}/cmr-matrix - -i sparse -o dense input.dense')
    call(f'{UNIMOD_DIRECTORY}/unimodularity-test modwheel_input.dense -s 2> /dev/null | egrep \'^[ 0-9-]*$\' 1> modwheel_signed.dense')
    call(f'{UNIMOD_DIRECTORY}/unimodularity-test modwheel_signed.dense -t -v 1> {file_base}-unimod.out 2> {file_base}-unimod.err')
    call(f'{UNIMOD_DIRECTORY}/unimodularity-test modwheel_signed.dense -t -v 1> {file_base}-unimodcert.out 2> {file_base}-unimodcert.err')

    # Remove temporary files.
    call(f'rm modwheel_input.dense  modwheel_signed.dense')

    r += 1

