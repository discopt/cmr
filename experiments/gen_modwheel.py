import sys
import os
import math
import subprocess
from config_modwheel import *

assert os.path.exists(INSTANCE_DIRECTORY)

for sample in SAMPLES:
  for instance in INSTANCES:
    file_base = f'{INSTANCE_DIRECTORY}/modwheel-{instance:05d}x{instance:05d}#{sample:03d}'
    os.system(f'{BUILD_DIRECTORY}/cmr-generate-wheel -01 {instance} -o sparse | {BUILD_DIRECTORY}/cmr-matrix -i sparse - -r -R2 {instance//2+1} -o sparse - | gzip > {file_base}.sparse.gz')

