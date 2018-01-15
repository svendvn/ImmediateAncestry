from numpy.random import choice
from proposal_permutation import Permutation
from proposal_substitution import Substitution
from proposal import Proposal

proposals=[Permutation, Substitution]


class MetaProposal(Proposal):
    
    def __init__(self, pops=[]):
        super(MetaProposal, self).__init__(pops)
        self.props=[]
        for prop_class in proposals:
            self.props.append(prop_class(pops=pops))
        
        
    def __call__(self, x, pks={}):
        i=choice(len(proposals), 1)
        prop=self.props[i]
        newx, g1,g2=prop(x)
        return newx, g1,g2,1,1,1