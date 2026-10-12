import os
import sys
import re
import graphlib
from local import *

TIME_LIMIT = 3600

runs = []

def print_usage(message=None):
  if message is not None:
    sys.stderr.write(f'{message}\n')
  sys.stderr.write(f'Usage: {sys.argv[0]} PATTERN [COMMAND]\n')
  sys.stderr.write('Shows or deletes runs matching PATTERN.\n')
  sys.stderr.write('Commands:\n')
  sys.stderr.write('  print   Prints the matching runs (default).\n')
  sys.stderr.write('  status  Prints the matching runs along with status information.\n')
  sys.stderr.write('  delete  Deletes the matching runs.\n')

def add_identity(num_samples):
  for sample in range(1, num_samples+1):
    for order in [1000, 2000, 4000]:
      instance = f'identity-{order:05}x{order:05}-{sample:03d}-gen'
      yield (instance, [])
      yield (f'identity-{order:05}x{order:05}-{sample:03d}-cmrdec', [instance])
      yield (f'identity-{order:05}x{order:05}-{sample:03d}-unimod', [instance])

def add_modwheel(num_samples):
  for sample in range(1, num_samples+1):
    for order in [13, 15, 17, 19, 21, 23, 25, 27, 101, 151, 201, 251, 301, 401, 501, 601, 701]:
      instance = f'modwheel-{order:05}x{order:05}-{sample:03d}-gen'
      yield (instance, [])
      if order <= 27:
        yield (f'modwheel-{order:05}x{order:05}-{sample:03d}-cmreuler', [instance])
        yield (f'modwheel-{order:05}x{order:05}-{sample:03d}-cmrpart', [instance])
      if order <= 251 or order == 301 and sample <= 5:
        yield (f'modwheel-{order:05}x{order:05}-{sample:03d}-cmrcert', [instance])
      if order <= 151 or order == 201 and sample <= 5:
        yield (f'modwheel-{order:05}x{order:05}-{sample:03d}-unimodcert', [instance])
      if order <= 701:
        yield (f'modwheel-{order:05}x{order:05}-{sample:03d}-cmrdec', [instance])
      if order <= 701:
        yield (f'modwheel-{order:05}x{order:05}-{sample:03d}-unimod', [instance])

def add_network(num_samples):
  for sample in range(1, num_samples+1):
    for order in [19, 20, 21, 22, 23, 24, 25, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 10000, 20000, 30000, 40000, 50000]:
      instance = f'network-{order:05}x{order:05}-{sample:03d}-gen'
      yield (instance, [])
      if order <= 30:
        yield (f'network-{order:05}x{order:05}-{sample:03d}-cmreuler', [instance])
        yield (f'network-{order:05}x{order:05}-{sample:03d}-cmrpart', [instance])
      if order <= 2500:
        yield (f'network-{order:05}x{order:05}-{sample:03d}-unimod', [instance])
      yield (f'network-{order:05}x{order:05}-{sample:03d}-cmrcamion', [instance])
      yield (f'network-{order:05}x{order:05}-{sample:03d}-cmrdec', [instance])

def add_rndcamion(num_samples):
  for sample in range(1, num_samples+1):
    for order in [100, 200, 300, 400, 500, 600, 700, 800]:
      for p in ['0.67']:
        instance = f'rndcamion-{order:05}x{order:05}-{p}-{sample:03d}-gen'
        yield (instance, [])
        yield (f'rndcamion-{order:05}x{order:05}-{p}-{sample:03d}-cmreuler', [instance])
        yield (f'rndcamion-{order:05}x{order:05}-{p}-{sample:03d}-cmrpart', [instance])
        yield (f'rndcamion-{order:05}x{order:05}-{p}-{sample:03d}-cmrdec', [instance])
        yield (f'rndcamion-{order:05}x{order:05}-{p}-{sample:03d}-unimod', [instance])
        yield (f'rndcamion-{order:05}x{order:05}-{p}-{sample:03d}-cmrcert', [instance])
        yield (f'rndcamion-{order:05}x{order:05}-{p}-{sample:03d}-unimodcert', [instance])


runs.extend(add_modwheel(100))
runs.extend(add_network(100))
runs.extend(add_rndcamion(100))
runs.extend(add_identity(1))

# Take care of dependencies.
runs_predecessors = { run[0]:run[1] for run in runs }
ts = graphlib.TopologicalSorter()
for run in runs:
  ts.add(run[0], *run[1])
runs = list(ts.static_order())

if __name__ == '__main__':

  if len(sys.argv) >= 2:
    re_run = re.compile(sys.argv[1])
  else:
    re_run = re.compile('.*')

  if len(sys.argv) >= 3:
    command = sys.argv[2]
  else:
    command = None

  for run in runs:
    re_result = re_run.match(run)
    if re_result:
      progress_file = f'{RESULT_STORAGE}/{run}.progress'
      if command == 'delete':
        if os.path.exists(progress_file):
          print(f'Deleting {run}')
          os.remove(progress_file)
      elif command == 'status':
        exists = os.path.exists(progress_file)
        if not exists:
          status = 'MISSING'
        else:
          progress_data = open(progress_file, 'r').read()
          status = progress_data.split('\n')[0]
        print(f'{run} [{status}]')
      elif command is None or command == 'print':
        print(run)
      else:
        print_usage(f'Invalid command <{command}>')
        sys.exit(1)


