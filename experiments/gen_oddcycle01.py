import sys
import os
import math
import subprocess

BUILD_DIRECTORY = '../build-release'
OUTPUT_DIRECTORY = 'oddcycle01'

try:
  num_repetitions = int(sys.argv[1])
except:
  print(f'Usage: {sys.argv[0]} NUM-REPETITIONS')
  sys.exit(1)

if not os.path.exists(OUTPUT_DIRECTORY):
  os.mkdir(OUTPUT_DIRECTORY)

for order in [11, 13, 15, 17, 19, 21, 23, 51, 101, 151, 201, 251, 301, 401, 501, 601, 701]:
  for r in range(1, num_repetitions+1):
    os.system(f'{BUILD_DIRECTORY}/cmr-generate-cycle -01 {order} -o sparse | {BUILD_DIRECTORY}/cmr-matrix -i sparse - -r -R2 {order//2} -o sparse {OUTPUT_DIRECTORY}/oddcycle01-{order:04d}x{order:04d}#{r:03d}.sparse')
