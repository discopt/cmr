import sys
import re
from math import isfinite

class Result:
  def __init__(self, instance, algorithm, sample=None):
    self._sample = sample
    self._instance = instance
    self._algorithm = algorithm
    self.time = float('inf')
    self._read = False

  @property
  def read(self):
    return self._read

  @property
  def instance(self):
    return self._instance

  @property
  def algorithm(self):
    return self._algorithm

  @property
  def sample(self):
    return self._sample

  def __repr__(self):
    sample_str = '' if self._sample is None else f' (sample #{self._sample})'
    return f'{self._instance} via {self._algorithm}{sample_str}: T={self.time}'

class BoxPlot:
  def __init__(self, instance, algorithm, data):
    self._instance = instance
    self._algorithm = algorithm
    n = len(data)
    self._num_finite = 0
    self._num_data = n
    for d in data:
      if isfinite(d):
        self._num_finite += 1
    self._average = sum(data) / n
    sorted_data = sorted(data)
    if len(sorted_data) % 2 == 1:
      self._median = sorted_data[n//2]
    else:
      self._median = 0.5 * sorted_data[n//2] + 0.5 * sorted_data[n//2-1]
    self._Q1 = sorted_data[n//4]
    self._Q3 = sorted_data[n - n//4 - 1]
    self._min = sorted_data[0]
    self._max = sorted_data[-1]

  @property
  def average(self):
    return self._average

  @property
  def median(self):
    return self._median

  @property
  def Q1(self):
    return self._Q1

  @property
  def Q3(self):
    return self._Q3

  @property
  def quantiles(self):
    return (self._Q1, self._median, self._Q3)

  @property
  def min(self):
    return self._min

  @property
  def max(self):
    return self._max

  @property
  def num_data(self):
    return self._num_data

  @property
  def num_finite(self):
    return self._num_finite

  def __repr__(self):
    return f'<bar plot: quantiles {self.Q1},{self.median},{self.Q3}, range [{self.min}, {self.max}], average {self.average}, {self.num_finite}/{self.num_data} finite>'

def average_times(result_list):
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

def compute_box_plots(result_list):
  result_map = { (result.instance, result.algorithm): [] for result in result_list }
  for result in result_list:
    if result.read:
      result_map[result.instance,result.algorithm].append(result.time)

  return { key: BoxPlot(key[0], key[1], times) for key,times in result_map.items() if times }

def parse_CMR(file_base, instance, algorithm, sample):
  result = Result(instance, algorithm, sample)

  sys.stderr.write(f'Reading {file_base}.out|.err\n')
  sys.stderr.flush()

  re_part_time = re.compile('  partition total: [0-9]* in ([0-9.]*) seconds')
  re_euler_time = re.compile('  eulerian enumeration total: [0-9]* in ([0-9.]*) seconds')
  re_dec_time = re.compile('  seymour total: [0-9]* in ([0-9.]*) seconds')
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
      result._read = True
  except:
    pass

  return result

def parse_network_camion(file_base, instance, algorithm, sample):
  result = Result(instance, algorithm, sample)

  sys.stderr.write(f'Reading {file_base}.out|.err\n')
  sys.stderr.flush()

  re_camion_time = re.compile('  general: [0-9]* in ([0-9.]*) seconds')
  re_graphic_time = re.compile('  total: [0-9]* in ([0-9.]*) seconds')
  re_time_limit = re.compile('Time limit exceeded!')

  time_limit = False
  camion_time = 0.0
  graphic_time = 0.0
  try:
    with open(file_base + '.err', 'r') as err_file:
      for line in err_file:
        match = re_time_limit.match(line)
        if match:
          time_limit = True
          continue 
        match = re_camion_time.match(line)
        if match and not time_limit:
          camion_time = float(match.group(1))
          continue
      result._read = True
    if result._read:
      result._read = False
      with open(file_base[:-9] + 'cmrgraphic.err', 'r') as err_file:
        for line in err_file:
          match = re_time_limit.match(line)
          if match:
            time_limit = True
            continue 
          match = re_graphic_time.match(line)
          if match and not time_limit:
            graphic_time = float(match.group(1))
            continue
        result._read = True
  except:
    pass

  result.time = camion_time + graphic_time

  return result

def parse_unimodularity_test(file_base, instance, algorithm, sample):
  result = Result(instance, algorithm, sample)

  sys.stderr.write(f'Reading {file_base}.out\n')

  re_total_time = re.compile('Total time: ([0-9.]*)')

  try:
    with open(file_base + '.out', 'r') as out_file:
      for line in out_file:
        match = re_total_time.match(line)
        if algorithm in ['unimod', 'unimodcert'] and match:
          result.time = float(match.group(1))
          continue
      result._read = True
  except:
    pass

  return result

def parse(file_base, instance, algorithm, sample):
  if algorithm[:9] == 'cmrcamion':
    return parse_network_camion(file_base, instance, algorithm, sample)
  elif algorithm[:3] == 'cmr':
    return parse_CMR(file_base, instance, algorithm, sample)
  elif algorithm[:6] == 'unimod':
    return parse_unimodularity_test(file_base, instance, algorithm, sample)

if __name__ == '__main__':
  result = parse(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))
  print(result)

