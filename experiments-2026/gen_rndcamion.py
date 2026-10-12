import sys
import os
import math
import subprocess

from config_rndcamion import *

assert os.path.exists(INSTANCE_DIRECTORY)

for sample in SAMPLES:
  for instance in INSTANCES:
    order,p = instance
    file_base = f'{INSTANCE_DIRECTORY}/rndcamion-{order:05d}x{order:05d}-p{p:0.2f}#{sample:03d}'
    os.system(f'{BUILD_DIRECTORY}/cmr-generate-random -o sparse {order} {order} {p} | {BUILD_DIRECTORY}/cmr-camion -i sparse - -S - | gzip > {file_base}.sparse.gz')

