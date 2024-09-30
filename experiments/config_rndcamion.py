from config import *

INSTANCE_DIRECTORY = 'rndcamion'
INSTANCES = [ (order,prob) for order in range(50,2001,50) for prob in [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1] ]
SAMPLES = list(range(1, 10+1))

