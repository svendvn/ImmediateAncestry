from numpy.random import choice
from copy import deepcopy
from proposal import Proposal

class Substitution(Proposal):
    
    def __call__(self, config):
        print("__call__(self, config). Substitution(Proposal)")

        config=list(config)
        print ("config", config)
        i=choice(len(config))
        print ("i", i)
        bases=deepcopy(self.pops)
        print ( "bases = ", bases)
        #print('config', config)
        #print(bases, config[i])
        bases.remove(config[i])
        print ( "bases (remove config[i])= ", bases)
        config2=deepcopy(config)
        #print(choice(bases,1))
        config2[i]=choice(bases,1)[0]
        #print(config2)
        return "".join(config2),1,1
    
if __name__=='__main__':
    s=Substitution(pops=list('evst'))
    s('evsvvstst')
