import sys
import os
import math
import subprocess

from config import *
OUTPUT_DIRECTORY = 'rndcamion'

try:
  num_repetitions = int(sys.argv[1])
except:
  print(f'Usage: {sys.argv[0]} NUM-REPETITIONS')
  sys.exit(1)

if not os.path.exists(OUTPUT_DIRECTORY):
  os.mkdir(OUTPUT_DIRECTORY)

for r in range(1, num_repetitions+1):
  for order in list(range(50, 2001, 50)):
    for p in [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]:
      file_base = f'{OUTPUT_DIRECTORY}/rndcamion-{order:05d}x{order:05d}-p{p:0.2f}#r{r:04d}'
      os.system(f'{BUILD_DIRECTORY}/cmr-generate-random -o sparse {order} {order} {p} | {BUILD_DIRECTORY}/cmr-camion -i sparse - -S - | gzip > {file_base}.sparse.gz')

