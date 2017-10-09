from numpy.random import choice
from copy import deepcopy

def permutation(config):
    if len(set(list(config)))==1:
        return config,1,1
    else:
        i,j=choice(range(len(config)),2,replace=False)
        a,b=config[i],config[j]
        while a==b:
            i,j=choice(range(len(config)),2,replace=False)
            a,b=config[i],config[j]
        config2=deepcopy(config)
        config2[i],config[j]=config[j],config[i]
    return config2, 1,1
    