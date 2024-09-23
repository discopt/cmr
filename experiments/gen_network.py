import sys
import os
import math
import subprocess

from config import *
OUTPUT_DIRECTORY = 'network'

try:
  num_repetitions = int(sys.argv[1])
except:
  print(f'Usage: {sys.argv[0]} NUM-REPETITIONS')
  sys.exit(1)

if not os.path.exists(OUTPUT_DIRECTORY):
  os.mkdir(OUTPUT_DIRECTORY)

for r in range(1, num_repetitions+1):
  for order in sorted(list(range(1,41)) + [ 100 * i for i in range(1,41) ] + [ 1000 * i for i in range(1,41) ]):
    file_base = f'{OUTPUT_DIRECTORY}/network-{order:05d}x{order:05d}#{r:03d}'
    os.system(f'{BUILD_DIRECTORY}/cmr-generate-network {order} {order} -o sparse | gzip > {file_base}.sparse.gz')

