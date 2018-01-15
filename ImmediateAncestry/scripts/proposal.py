class Proposal(object):
    
    def __init__(self, pops):
        self.pops=pops
    
    def __call__(self,x, pks):
        return x,1,1,1,1,1
    
        