'''
Created on 11/10/2016

@author: svendvn
'''


import tmhmm
import numpy
import matrix_generation
from math import log

import operator
import functools

from id_dic import id_dic

from itertools import product

import warnings

code_to_index={"A":0, "C":1, "G":2, "T":3}


# def call_helper():
#     print(help(tmhmm.hmm.forward))
#     
# def call_forward():
#     init=numpy.array([.5,.5])
#     transes=numpy.array([[[.5,.5],[.2,.8]], [[.1,.9],[.4,.6]]])
#     emiss=numpy.array([[.9,.1],[.95,.05]])
#     sequence="".join(["a","b","b","b","b","a","a","a","b","a"])
#     char_map={"A":0, "B":1}
#     def to_call_each(i):
#         return 0
#     
#     return tmhmm.hmm.forward(sequence, 
#                              to_call_each,
#                              init,
#                              transes,
#                              emiss,
#                              char_map,
#                              None,
#                              None) 
    
class likelihood_class(object):
    
    def __init__(self, transition_generator, emission_generator_function, initial_probabilities, sequence, char_map, short_to_full):
        self.transition_generator=transition_generator
        self.emission_generator_function=emission_generator_function
        self.initial_probabilities=initial_probabilities
        self.sequence=sequence
        self.char_map=char_map
        self.short_to_full=short_to_full
        self.values_dic={}
        
    def __call__(self, params, pks={}):
        
        reduced=tuple(find_smallet_equivalence_class(params))
        if reduced in self.values_dic:
            return self.values_dic[reduced]
        else:
            params2=[self.short_to_full[p] for p in params]
            emission_generator=self.emission_generator_function(params2)
            _, Cs= tmhmm.hmm.forward2(self.sequence, numpy.array(self.initial_probabilities), self.transition_generator, emission_generator, self.char_map, None, None)
            print(Cs)
            res=0
            for C in Cs:
                if C<=0:
                    self.values_dic[reduced]=-float('Inf')
                    return -float('Inf')
                res=res+log(C)
            self.values_dic[reduced]=res
            return res
    
# def generate_likelihood_from_generators(transition_generator, emission_generator_function, initial_probabilities, sequence, char_map):
#     def likelihood(params):
#         #return 0 #for testing
#         emission_generator=emission_generator_function(params)
#         
#         _, Cs= tmhmm.hmm.forward2(sequence, numpy.array(initial_probabilities), transition_generator, emission_generator, char_map, None, None)
#         ##FIXME: math domain error
#         res=0
#         for C in Cs:
#             if C<=0:
#                 return -float('Inf')
#             res=res+log(C)
#         return res
#     return likelihood



def generate_binned_likelihood_from_data(bin_maps, alleles_list, recomb_map_list, seq_list, possible_pops, generations=3, short_to_full=id_dic, rho_infinity=False):
    list_of_likelihoods=[]
    for bin_map, alleles, recomb_map, seq in zip(bin_maps, alleles_list, recomb_map_list, seq_list):
        trans_gen=matrix_generation.generate_transition_matrix(recomb_map, generations, rho_infinity)
        before_dic=matrix_generation.calculate_before_dic(bin_map, alleles, seq, possible_pops)
        ems_gen=matrix_generation.generate_emission_matrix_binned(before_dic, generations)
        M=(2**(generations-1))**2
        initial=[1.0/M]*M
        char_map={str(i):i for i in range(6)}
        seq2="0"*len(bin_map)
        if len(bin_map)<2:
            message='length of bin_map is just '+str(len(bin_map))
            warnings.warn(message, UserWarning)
        list_of_likelihoods.append(likelihood_class(trans_gen, ems_gen, initial, seq2, char_map, short_to_full))
        
    def likelihood(param,pks={}):
        return sum([lik(param,pks) for lik in list_of_likelihoods])
    
    return likelihood
    

def generate_likelihood_from_data(alleles_list, recomb_map_list, seq_list, generations=3, short_to_full=id_dic(), rho_infinity=False):
    list_of_likelihoods=[]
    for alleles, recomb_map, seq in zip(alleles_list, recomb_map_list, seq_list):
        trans_gen=matrix_generation.generate_transition_matrix(recomb_map, generations, rho_infinity)
        ems_gen=matrix_generation.generate_emission_matrix(alleles, generations)
        M=(2**(generations-1))**2
        initial=[1.0/M]*M
        char_map={str(i):i for i in range(6)}
        seq2="".join(map(str,seq))
        list_of_likelihoods.append(likelihood_class(trans_gen, ems_gen, initial, seq2, char_map, short_to_full))
        
    def likelihood(param,pks={}):
        return sum([lik(param,pks) for lik in list_of_likelihoods])
    
    return likelihood
    


# def test_model_likelihood(size):
#     alleles={"pop1":[0.2]*size, "pop2":[0.8]*size}
#     popsm=["pop2","pop2","pop2","pop2"]
#     popsf=["pop2","pop2","pop2","pop2"]
#     trans_gen=matrix_generation.generate_transition_matrix([0.1]*size)
#     mu=0.01
#     ems_gen=matrix_generation.generate_emission_matrix(alleles)
#     initial=[1/16.0]*16
#     major_allele="C"*size
#     #seq="CCMWCRSA"[:size]+"C"*(size-8)
#     seq=major_allele
#   
#     seq={key:transform_from_fasta_to_index(key,ma) for key,ma in zip(seq,major_allele)}
#     lik=generate_likelihood(trans_gen, ems_gen, initial, seq, char_map)
#     return lik
    
def transform_from_fasta_to_index(x, ma="C"):
    def transform(x):
        if x=='Y':
            return 1,0#"T","C"
        if x=='M':
            return 1,0#"A","C"
        if x=='S':
            return 0,1#"C","G"
        if x=="C":
            return 0,0
        return 1,1
    def transform2(tup):
        index=2*tup[0]+tup[1]
        return index
    return transform2(transform(x))

def get_val(listi):
    return "".join(listi)

def get_best_of_size(listi, size):
    assert len(listi)%size==0, "wrong list size"
    assert size%2==0, "wrong size"
    band=listi
    h= size//2
    bval= get_val(listi)
    for i in range(len(listi)//size):
        cand=band[:(size*i)]+\
                band[(size*i+h):(size*(i+1))]+\
                band[(size*i):(size*i+h)]+\
                band[(size*(i+1)):]
        
        cval=get_val(cand)
        if cval<bval:
            band=cand
            bval=cval
    return bval



def find_smallet_equivalence_class(params):
    size=2
    best_combi=params
    while size<=len(params):
        best_combi=get_best_of_size(best_combi, size)
        size*=2
    return best_combi
    
        


if __name__ == '__main__':
    #call_helper()
    #print(call_forward())
    #lik=test_model_likelihood(8)
    #print(maximize_likelihood_exhaustive(lik,["pop1","pop2"]))
    print(find_smallet_equivalence_class(["a","b","a","a","a","b","a","a"]))
    n=2**4
    combinations=[]
    counter=0
    for a,itera in enumerate(product(*([['a','b','c','d']]*n))):
        #print(itera)
        if a%10000==0:
            print(itera)
        if find_smallet_equivalence_class(itera) in combinations:
            #print(str(itera), "skipped")
            continue
        counter+=1
        #combinations.append(itera)
        
        #likelihoods.append(counter)
        #print(itera)
    print('counter',counter)
        
    
    