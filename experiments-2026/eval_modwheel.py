import sys
import os
import math
import subprocess
import eval_parse
import re
from local import *
from runs import runs, TIME_LIMIT

assert os.path.exists(RESULT_STORAGE)

run_filter = re.compile('.*' if len(sys.argv) <= 1 else sys.argv[1])
command = None if len(sys.argv) <= 2 else sys.argv[2]

results = []
for run in runs:
  name = run.split('-')
  if name[0] != 'modwheel' or name[-1] == 'gen':
    continue

  match_filter = run_filter.search(name[1])
  if not match_filter:
    continue

  result = eval_parse.parse(RESULT_STORAGE + '-'.join(name), '-'.join(name[:-2]), name[-1], int(name[-2]))
  if result is None:
    sys.stderr.write(f'No result for {run}.\n')
    sys.stderr.flush()
    sys.exit(1)

  results.append(result)
  print(f'    {result}')

sys.stderr.write(f'Parsed {len(results)} logs.\n')
sys.stderr.flush()

times = eval_parse.compute_box_plots(results)

sys.stderr.write(f'Extracted {len(times)} entries.\n')
sys.stderr.flush()

if command is not None and command.startswith('boxplot'):
  # `boxplot` command.
  for key,plot in times.items():
    print(f'% {key[1]} for {key[0]} with {plot.num_data} samples')
    def limit_time(t):
      if t >= TIME_LIMIT:
        return 2*TIME_LIMIT
      else:
        return t
    print(f'\\addplot[boxplot-{key[1]}, boxplot prepared={{draw position=TODO,%')
    print(f'    median={limit_time(plot.median)}, upper quartile={limit_time(plot.Q3)}, lower quartile={limit_time(plot.Q1)}, upper whisker={limit_time(plot.max)}, lower whisker={limit_time(plot.min)}}}] coordinates {{}};')
else:
  # Others: just print.
  for key,plot in times.items():
    print(f'{key[0]} via {key[1]}: {plot}')

