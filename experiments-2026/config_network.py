from config import *

INSTANCE_DIRECTORY = 'network'
INSTANCES = sorted(set(range(1,41)) | set(100 * i for i in range(1,41) ) | set(1000 * i for i in range(1,41) ))
SAMPLES = list(range(1, 10+1))

