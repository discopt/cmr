import sys
import os
import math
import subprocess
import eval_parse

from config_modwheel import *

assert os.path.exists(INSTANCE_DIRECTORY)

results = []
for instance in INSTANCES:
  file_prefix = f'{INSTANCE_DIRECTORY}/modwheel-{instance:05d}x{instance:05d}'
  for algo in ALGORITHMS:
    for sample in SAMPLES:
      result = eval_parse.parse(file_prefix, instance, algo, sample)

      if algo.startswith('unimod') and result.time >= 3700.0:
        result.time = float('inf')

      results.append( result )
      print(f'    {result}')

times = eval_parse.averageTimes(results)

for instance in INSTANCES:
  for algo in ALGORITHMS:
    print(instance, algo, times[ (instance, algo) ])

