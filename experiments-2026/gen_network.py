import sys
import os
import math
import subprocess

from config_network import *

assert os.path.exists(INSTANCE_DIRECTORY)

for sample in SAMPLES:
  for instance in INSTANCES:
    order = instance
    file_base = f'{INSTANCE_DIRECTORY}/network-{order:05d}x{order:05d}#{sample:03d}'
    os.system(f'{BUILD_DIRECTORY}/cmr-generate-network {order} {order} -o sparse | gzip > {file_base}.sparse.gz')

