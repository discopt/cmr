import sys
import re

class Result:
  def __init__(self, instance, algorithm, run=None):
    self._run = run
    self._instance = instance
    self._algorithm = algorithm
    self.time = float('inf')
    self.time_max = 0.0
    self.time_min = 0.0

  @property
  def instance(self):
    return self._instance

  @property
  def algorithm(self):
    return self._algorithm

  @property
  def run(self):
    return self._run

  def __repr__(self):
    run_str = '' if self._run is None else f' (run #{self._run})'
    return f'{self._instance} via {self._algorithm}{run_str}: {self.time} in [{self.time_min},{self.time_max}]'

def averageTimes(result_list):
  result_map = { (result.instance, result.algorithm): [] for result in result_list }
  for result in result_list:
    result_map[result.instance,result.algorithm].append(result)

  def avg(instance, algorithm, L):
    times = [ result.time for result in L ]
    res = Result(instance, algorithm)
    res.time = sum(times) / len(L)
    res.time_min = min(times)
    res.time_max = max(times)
    return res

  return { key: avg(key[0], key[1], value) for key,value in result_map.items() }

def parseCMR(file_base, instance, algorithm, run):
  result = Result(instance, algorithm, run)

  sys.stderr.write(f'Reading {file_base}.out|.err\n')

  re_part_time = re.compile('  partition time: ([0-9.]*)( seconds|)')
  re_euler_time = re.compile('  eulerian enumeration time: ([0-9.]*)( seconds|)')
  re_dec_time = re.compile('  seymour total: [0-9]* in ([0-9.]*)( seconds|)')
  re_time_limit = re.compile('Time limit exceeded!')

  try:
    time_limit = False
    with open(file_base + '.err', 'r') as err_file:
      for line in err_file:
        match = re_time_limit.match(line)
        if match:
          time_limit = True
          continue 
        match = re_part_time.match(line)
        if algorithm == 'cmrpart' and match and not time_limit:
          result.time = float(match.group(1))
          continue
        match = re_euler_time.match(line)
        if algorithm == 'cmreuler' and match and not time_limit:
          result.time = float(match.group(1))
          continue
        match = re_dec_time.match(line)
        if algorithm in ['cmrdec', 'cmrcert'] and match and not time_limit:
          result.time = float(match.group(1))
          continue
  except:
    pass

  return result

def parseUnimodularityTest(file_base, instance, algorithm, run):
  result = Result(instance, algorithm, run)

  sys.stderr.write(f'Reading {file_base}.out\n')

  re_total_time = re.compile('Total time: ([0-9.]*)')

  try:
    with open(file_base + '.out', 'r') as out_file:
      for line in out_file:
        match = re_total_time.match(line)
        if algorithm in ['unimod', 'unimodcert'] and match:
          result.time = float(match.group(1))
          continue
  except:
    pass

  return result

def parse(file_prefix, instance, algorithm, run):
  if algorithm[:3] == 'cmr':
    return parseCMR(f'{file_prefix}#{run:03d}-{algorithm}', instance, algorithm, run)
  elif algorithm[:6] == 'unimod':
    return parseUnimodularityTest(f'{file_prefix}#{run:03d}-{algorithm}', instance, algorithm, run)

if __name__ == '__main__':
  result = parse(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))
  print(result)

