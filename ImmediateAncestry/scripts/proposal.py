class Proposal(object):
    
    def __init__(self, pops):
        print ("#__init__(self, pops)   Proposal(object)")
        self.pops=pops
        print ("\n")
    
    def __call__(self,x, pks):
        print ("#__call__(self, x,pks)")
        print ("return x,,1,1,1,1,1 = ", x,1,1,1,1,1 , "\n")
        return x,1,1,1,1,1
    
        
