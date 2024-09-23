import sys
import os
import math
import subprocess
import eval_parse

from config import *
OUTPUT_DIRECTORY = 'modwheel'

try:
  num_repetitions = int(sys.argv[1])
except:
  print(f'Usage: {sys.argv[0]} NUM-REPETITIONS')
  sys.exit(1)

if not os.path.exists(OUTPUT_DIRECTORY):
  os.mkdir(OUTPUT_DIRECTORY)

results = []
for order in sorted(list(range(9, 100, 10)) + list(range(99, 1000, 100))):
  print(f'{order}')
  file_prefix = f'{OUTPUT_DIRECTORY}/modwheel-{order:05d}x{order:05d}'
  for algo in ['cmrdec', 'cmrpart', 'cmreuler', 'cmrcert', 'unimod', 'unimodcert']:
    print(f'  {algo}')
    for r in range(1, num_repetitions+1):
      result = eval_parse.parse(file_prefix, order, algo, r)
      results.append( result )
      print(f'    {result}')

times = eval_parse.averageTimes(results)

for order in sorted(list(range(9, 99, 10)) + list(range(99, 1000, 100))):
  for algo in ['cmrdec', 'cmrpart', 'cmreuler', 'cmrcert', 'unimod', 'unimodcert']:
    print(order, algo, times[ (order, algo) ])


