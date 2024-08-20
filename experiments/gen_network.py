import sys
import os
import math
import subprocess

BUILD_DIRECTORY = '../build-release'
OUTPUT_DIRECTORY = 'network'

try:
  num_repetitions = int(sys.argv[1])
except:
  print(f'Usage: {sys.argv[0]} NUM-REPETITIONS')
  sys.exit(1)

if not os.path.exists(OUTPUT_DIRECTORY):
  os.mkdir(OUTPUT_DIRECTORY)

for order in [10, 12, 14, 16, 18, 20, 22, 24, 50, 100, 200, 400, 800, 1600]:
  for r in range(1, num_repetitions+1):
    output_file = open(f'{OUTPUT_DIRECTORY}/network-{order:04d}x{order:04d}#{r:03d}.sparse', 'w')
    process = subprocess.run([f'{BUILD_DIRECTORY}/cmr-generate-network', '-o', 'sparse', str(order), str(order)], stdout=output_file)
    output_file.close()
