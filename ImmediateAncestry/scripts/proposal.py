from numpy.random import choice
from proposal_permutation import permutation
from proposal_substitution import substitution

proposals=[permutation, substitution]

def proposal(x):
    i=choice(len(proposals), 1)
    prop=proposals[i]
    newx, g1,g2=prop(x)
    return newx, g1,g2,1,1,1