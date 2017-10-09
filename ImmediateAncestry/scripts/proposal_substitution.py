from numpy.random import choice
from copy import deepcopy

def substitution(config):
    i=choice(len(config))
    bases=['A','C','G','T']
    bases.remove(config[i])
    config2=deepcopy(config)
    config2[i]=choice(bases,1)
    return config2,1,1