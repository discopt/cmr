import re
import sys
import os
import time
from local import *
from runs import runs, runs_predecessors, TIME_LIMIT
from datetime import datetime
from filelock import Timeout, FileLock

pid = os.getpid()
temp_storage = f'{TEMP_STORAGE}/{pid}'
os.system(f'mkdir -p {temp_storage}')
  
# Argument 1: Filter for jobs.
if len(sys.argv) >= 2:
  re_run = re.compile(sys.argv[1])
else:
  re_run = re.compile('.*')

def call(command):
#  print(f'`{command}`')
  os.system(command)
      
print(f'Execution {pid}: starting loop at {datetime.now()}.', flush=True)

for run in runs:
  re_result = re_run.match(run)
  if not re_result:
    continue

  # Ensure we're the only process looking at this run.
  lock = FileLock(f'{RESULT_STORAGE}/{run}.lock', timeout=0)
  try:
    lock.acquire()
    with lock:

      # Is that run already done?
      if os.path.exists(f'{RESULT_STORAGE}/{run}.progress'):
        continue

      # Are all dependencies already done? If not, some other job must be processing a dependency and will look at this one, too.
      for pred in runs_predecessors[run]:
        while not os.path.exists(f'{RESULT_STORAGE}/{pred}.progress'):
          print(f'Execution {pid}: postponing {run} at {datetime.now()} to wait for {pred}.', flush=True)
          time.sleep(60)
          continue

      # Mark the run as running in the progress file.
      with open(f'{RESULT_STORAGE}/{run}.progress', 'w') as f:
        f.write('RUNNING\n')
        
      # We're good to go!
      print(f'Execution {pid}: running {run} at {datetime.now()}.', flush=True)
      name = run.split('-')
      
      if name[0] == 'identity':
        order = int(name[1].split('x')[0])
        sample = int(name[2])
        file_base = RESULT_STORAGE + '/' + '-'.join(name[:-1])
        instance_file = f'{file_base}.sparse.gz'
    
        if name[-1] == 'gen':
          call(f'{BUILD_DIRECTORY}/cmr-named -g I_{order} -o sparse - | gzip > {file_base}.sparse.gz')
    
        elif name[-1] == 'cmrdec':
          call(f'gunzip -cd {instance_file} | {BUILD_DIRECTORY}/cmr-tu - -i sparse --stats --algo decomposition --time-limit {TIME_LIMIT} 1> {file_base}-cmrdec.out 2> {file_base}-cmrdec.err')
    
        elif name[-1] == 'unimod':
          call(f'gunzip -cd {instance_file} | {BUILD_DIRECTORY}/cmr-matrix - -i sparse -o dense {temp_storage}/input.dense')
          call(f'{UNIMOD_DIRECTORY}/unimodularity-test {temp_storage}/input.dense -s 2> /dev/null | egrep \'^[ 0-9-]*$\' 1> {temp_storage}/signed.dense')
          call(f'{UNIMOD_DIRECTORY}/unimodularity-test {temp_storage}/signed.dense -t -v 1> {file_base}-unimod.out 2> {file_base}-unimod.err')
    
      if name[0] == 'modwheel':
    
        order = int(name[1].split('x')[0])
        sample = int(name[2])
        file_base = RESULT_STORAGE + '/' + '-'.join(name[:-1])
        instance_file = f'{file_base}.sparse.gz'
    
        if name[-1] == 'gen':
          call(f'{BUILD_DIRECTORY}/cmr-generate-wheel -01 {order} -o sparse | {BUILD_DIRECTORY}/cmr-matrix -i sparse - -r -R3 {order//2+1} -o sparse - | gzip > {instance_file}')
    
        elif name[-1] == 'cmreuler':
          call(f'gunzip -cd {instance_file} | {BUILD_DIRECTORY}/cmr-tu - -i sparse --stats --algo eulerian --time-limit {TIME_LIMIT} 1> {file_base}-cmreuler.out 2> {file_base}-cmreuler.err')
    
        elif name[-1] == 'cmrpart':
          call(f'gunzip -cd {instance_file} | {BUILD_DIRECTORY}/cmr-tu - -i sparse --stats --algo partition --time-limit {TIME_LIMIT} 1> {file_base}-cmrpart.out 2> {file_base}-cmrpart.err')
    
        elif name[-1] == 'cmrdec':
          call(f'gunzip -cd {instance_file} | {BUILD_DIRECTORY}/cmr-tu - -i sparse --stats --algo decomposition --time-limit {TIME_LIMIT} 1> {file_base}-cmrdec.out 2> {file_base}-cmrdec.err')
    
        elif name[-1] == 'cmrcert':
          call(f'gunzip -cd {instance_file} | {BUILD_DIRECTORY}/cmr-tu - -i sparse --stats --algo decomposition --time-limit {TIME_LIMIT} -N {file_base}-cmrcert.sub 1> {file_base}-cmrcert.out 2> {file_base}-cmrcert.err')
    
        elif name[-1] == 'unimod':
          call(f'gunzip -cd {instance_file} | {BUILD_DIRECTORY}/cmr-matrix - -i sparse -o dense {temp_storage}/input.dense')
          call(f'{UNIMOD_DIRECTORY}/unimodularity-test {temp_storage}/input.dense -s 2> /dev/null | egrep \'^[ 0-9-]*$\' 1> {temp_storage}/signed.dense')
          call(f'{UNIMOD_DIRECTORY}/unimodularity-test {temp_storage}/signed.dense -t -v 1> {file_base}-unimod.out 2> {file_base}-unimod.err')
    
        elif name[-1] == 'unimodcert':
          call(f'gunzip -cd {instance_file} | {BUILD_DIRECTORY}/cmr-matrix - -i sparse -o dense {temp_storage}/input.dense')
          call(f'{UNIMOD_DIRECTORY}/unimodularity-test {temp_storage}/input.dense -s 2> /dev/null | egrep \'^[ 0-9-]*$\' 1> {temp_storage}/signed.dense')
          call(f'{UNIMOD_DIRECTORY}/unimodularity-test {temp_storage}/signed.dense -t -v -c 1> {file_base}-unimodcert.out 2> {file_base}-unimodcert.err')
    
      if name[0] == 'network':
        order = int(name[1].split('x')[0])
        sample = int(name[2])
        file_base = RESULT_STORAGE + '/' + '-'.join(name[:-1])
        instance_file = f'{file_base}.sparse.gz'
    
        if name[-1] == 'gen':
          call(f'{BUILD_DIRECTORY}/cmr-generate-network {order} {order} -o sparse | gzip > {file_base}.sparse.gz')
    
        elif name[-1] == 'cmreuler':
          call(f'gunzip -cd {instance_file} | {BUILD_DIRECTORY}/cmr-tu - -i sparse --stats --algo eulerian --time-limit {TIME_LIMIT} 1> {file_base}-cmreuler.out 2> {file_base}-cmreuler.err')
    
        elif name[-1] == 'cmrpart':
          call(f'gunzip -cd {instance_file} | {BUILD_DIRECTORY}/cmr-tu - -i sparse --stats --algo partition --time-limit {TIME_LIMIT} 1> {file_base}-cmrpart.out 2> {file_base}-cmrpart.err')
    
        elif name[-1] == 'cmrdec':
          call(f'gunzip -cd {instance_file} | {BUILD_DIRECTORY}/cmr-tu - -i sparse --stats --algo decomposition --time-limit {TIME_LIMIT} 1> {file_base}-cmrdec.out 2> {file_base}-cmrdec.err')
    
        elif name[-1] == 'cmrcamion':
          call(f'gunzip -cd {file_base}.sparse.gz | {BUILD_DIRECTORY}/cmr-matrix - -i sparse -o sparse - -c | {BUILD_DIRECTORY}/cmr-graphic - -i sparse --stats --time-limit {TIME_LIMIT} 1> {file_base}-cmrgraphic.out 2> {file_base}-cmrgraphic.err')
          call(f'gunzip -cd {file_base}.sparse.gz | {BUILD_DIRECTORY}/cmr-camion - -i sparse --stats --time-limit {TIME_LIMIT} 1> {file_base}-cmrcamion.out 2> {file_base}-cmrcamion.err')
    
        elif name[-1] == 'unimod':
          call(f'gunzip -cd {instance_file} | {BUILD_DIRECTORY}/cmr-matrix - -i sparse -o dense {temp_storage}/input.dense')
          call(f'{UNIMOD_DIRECTORY}/unimodularity-test {temp_storage}/input.dense -s 2> /dev/null | egrep \'^[ 0-9-]*$\' 1> {temp_storage}/signed.dense')
          call(f'{UNIMOD_DIRECTORY}/unimodularity-test {temp_storage}/signed.dense -t -v 1> {file_base}-unimod.out 2> {file_base}-unimod.err')
    
      if name[0] == 'rndcamion':
        order = int(name[1].split('x')[0])
        p = name[2]
        sample = int(name[3])
        file_base = RESULT_STORAGE + '/' + '-'.join(name[:-1])
        instance_file = f'{file_base}.sparse.gz'
    
        if name[-1] == 'gen':
          call(f'{BUILD_DIRECTORY}/cmr-generate-random -o sparse {order} {order} {p} | {BUILD_DIRECTORY}/cmr-camion -i sparse - -S - | gzip > {file_base}.sparse.gz')
    
        elif name[-1] == 'cmreuler':
          call(f'gunzip -cd {instance_file} | {BUILD_DIRECTORY}/cmr-tu - -i sparse --stats --algo eulerian --time-limit {TIME_LIMIT} 1> {file_base}-cmreuler.out 2> {file_base}-cmreuler.err')
    
        elif name[-1] == 'cmrpart':
          call(f'gunzip -cd {instance_file} | {BUILD_DIRECTORY}/cmr-tu - -i sparse --stats --algo partition --time-limit {TIME_LIMIT} 1> {file_base}-cmrpart.out 2> {file_base}-cmrpart.err')
    
        elif name[-1] == 'cmrdec':
          call(f'gunzip -cd {instance_file} | {BUILD_DIRECTORY}/cmr-tu - -i sparse --stats --algo decomposition --time-limit {TIME_LIMIT} 1> {file_base}-cmrdec.out 2> {file_base}-cmrdec.err')
    
        elif name[-1] == 'cmrcert':
          call(f'gunzip -cd {instance_file} | {BUILD_DIRECTORY}/cmr-tu - -i sparse --stats --algo decomposition --time-limit {TIME_LIMIT} -N {file_base}-cmrcert.sub 1> {file_base}-cmrcert.out 2> {file_base}-cmrcert.err')
    
        elif name[-1] == 'unimod':
          call(f'gunzip -cd {instance_file} | {BUILD_DIRECTORY}/cmr-matrix - -i sparse -o dense {temp_storage}/input.dense')
          call(f'{UNIMOD_DIRECTORY}/unimodularity-test {temp_storage}/input.dense -s 2> /dev/null | egrep \'^[ 0-9-]*$\' 1> {temp_storage}/signed.dense')
          call(f'{UNIMOD_DIRECTORY}/unimodularity-test {temp_storage}/signed.dense -t -v 1> {file_base}-unimod.out 2> {file_base}-unimod.err')
    
        elif name[-1] == 'unimodcert':
          call(f'gunzip -cd {instance_file} | {BUILD_DIRECTORY}/cmr-matrix - -i sparse -o dense {temp_storage}/input.dense')
          call(f'{UNIMOD_DIRECTORY}/unimodularity-test {temp_storage}/input.dense -s 2> /dev/null | egrep \'^[ 0-9-]*$\' 1> {temp_storage}/signed.dense')
          call(f'{UNIMOD_DIRECTORY}/unimodularity-test {temp_storage}/signed.dense -t -v -c 1> {file_base}-unimodcert.out 2> {file_base}-unimodcert.err')
      
      # Mark the run as done in the progress file.
      with open(f'{RESULT_STORAGE}/{run}.progress', 'w') as f:
        f.write('DONE\n')
  
  except Timeout:
    print(f'Execution {pid}: skipping {run} due to lock.', flush=True)

print(f'Execution {pid}: stopping loop at {datetime.now()}.', flush=True)
os.system(f'rm -r {temp_storage}')

