from ggrandparents_model import find_smallet_equivalence_class
from math import factorial
from itertools import combinations, combinations_with_replacement
import sys, math
from time import time, localtime, strftime
from itertools import product
import re 
import os
import subprocess
import numpy

def prior (x):

    small = find_smallet_equivalence_class(x) 
    string=list(small)
    var="['"+"', '".join(string)+"']"
    i=0

    for line in open("archivo2.txt"):
       if var in line:
  
           prior_num = line.split("	")[1]
          
    prior_prob=numpy.log(1/(int(prior_num)))
    print (prior_prob)
    return prior_prob




def posterior(likelihood, prior):

    def post(x,pks={}):
        return likelihood(x)+prior(x)

    return post




