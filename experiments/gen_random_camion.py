import sys
import os
import math
import subprocess

BUILD_DIRECTORY = '../build-release'
OUTPUT_DIRECTORY = 'random-camion'

try:
  num_repetitions = int(sys.argv[1])
except:
  print(f'Usage: {sys.argv[0]} NUM-REPETITIONS')
  sys.exit(1)

if not os.path.exists(OUTPUT_DIRECTORY):
  os.mkdir(OUTPUT_DIRECTORY)

for order in [200, 400, 800, 1600]:
  for p in [0.6667, 0.5, 0.25, 0.125]:
    for r in range(1, num_repetitions+1):
      os.system(f'{BUILD_DIRECTORY}/cmr-generate-random -o sparse {order} {order} {p} | {BUILD_DIRECTORY}/cmr-camion -i sparse - -S {OUTPUT_DIRECTORY}/random-camion-{order:04d}x{order:04d}-p{p:0.3f}#{r:03d}.sparse')
