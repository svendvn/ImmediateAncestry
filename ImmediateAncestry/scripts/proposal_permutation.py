from numpy.random import choice
from copy import deepcopy
from proposal import Proposal

class Permutation(Proposal):
    
    def __call__(self, config):
        print ( "#__call__(self, config) . Permutation(Proposal)")
        config=list(config)
        print ("config = ", config)
        if len(set(list(config)))==1:
            print (" if len(set(list(config)))==1")
            print ( "return config,1,1", config, "\n")
            return ''.join(config),1,1

        else:
            print (" if len(set(list(config))) |=1")
            i,j=choice(range(len(config)),2,replace=False)
            a,b=config[i],config[j]
            print ("i,j,a,b = ", i,j,a,b)

            while a==b:
                print ("while a=b")
                i,j=choice(range(len(config)),2,replace=False)
                a,b=config[i],config[j]
                print ("i,j,a,b = ", i,j,a,b)

            config2=deepcopy(config)
            config2[i],config[j]=config[j],config[i]
            print("return config2,1,1", config2,"\n")

        return ''.join(config2), 1,1
    
