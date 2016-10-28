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

code_to_index={"A":0, "C":1, "G":2, "T":3}


def call_helper():
    print(help(tmhmm.hmm.forward))
    
def call_forward():
    init=numpy.array([.5,.5])
    transes=numpy.array([[[.5,.5],[.2,.8]], [[.1,.9],[.4,.6]]])
    emiss=numpy.array([[.9,.1],[.95,.05]])
    sequence="".join(["a","b","b","b","b","a","a","a","b","a"])
    char_map={"A":0, "B":1}
    def to_call_each(i):
        return 0
    
    return tmhmm.hmm.forward(sequence, 
                             to_call_each,
                             init,
                             transes,
                             emiss,
                             char_map,
                             None,
                             None) 
    
def generate_likelihood_from_generators(transition_generator, emission_generator_function, initial_probabilities, sequence, char_map):
    def likelihood(params):
        emission_generator=emission_generator_function(params)
        _, Cs= tmhmm.hmm.forward2(sequence, numpy.array(initial_probabilities), transition_generator, emission_generator, char_map, None, None)
        return sum([log(C) for C in Cs])
    return likelihood

def maximize_likelihood_exhaustive(likelihood, pops_to_choose_from):
    combinations=[]
    likelihoods=[]
    counter=0
    for fff in pops_to_choose_from:
        for ffm in pops_to_choose_from:
            for fmf in pops_to_choose_from:
                for fmm in pops_to_choose_from:
                    for mff in pops_to_choose_from:
                        for mfm in pops_to_choose_from:
                            for mmf in pops_to_choose_from:
                                for mmm in pops_to_choose_from:
                                    itera=(fff,ffm,fmf,fmm,mff,mfm,mmf,mmm)
                                    if find_smallet_equivalence_class(itera) in combinations:
                                        print(str(itera), "skipped")
                                        continue
                                    counter+=1
                                    combinations.append(itera)
                                    #likelihoods.append(counter)
                                    likelihoods.append(likelihood(list(itera)))
    sorted_indexes=[i[0] for i in sorted(enumerate(likelihoods), key=lambda x: x[1])]
    res_dict=[]
    for i in range(10):
        index=sorted_indexes[-i-1]
        res_dict.append((combinations[index],likelihoods[index]))
    print("Looped over ",counter)
    return res_dict
    

def generate_likelihood_from_data(alleles_list, recomb_map_list, seq_list):
    list_of_likelihoods=[]
    for alleles, recomb_map, seq in zip(alleles_list, recomb_map_list, seq_list):
        trans_gen=matrix_generation.generate_transition_matrix(recomb_map)
        ems_gen=matrix_generation.generate_emission_matrix(alleles)
        initial=[1/16.0]*16
        char_map={'0':0,'1':1,'2':2}
        seq2="".join(map(str,seq))
        list_of_likelihoods.append(generate_likelihood_from_generators(trans_gen, ems_gen, initial, seq2, char_map))
        
    def likelihood(param):
        return sum([lik(param) for lik in list_of_likelihoods])
    
    return likelihood
    


def test_model_likelihood(size):
    alleles={"pop1":[0.2]*size, "pop2":[0.8]*size}
    popsm=["pop2","pop2","pop2","pop2"]
    popsf=["pop2","pop2","pop2","pop2"]
    trans_gen=matrix_generation.generate_transition_matrix([0.1]*size)
    mu=0.01
    ems_gen=matrix_generation.generate_emission_matrix(alleles)
    initial=[1/16.0]*16
    major_allele="C"*size
    #seq="CCMWCRSA"[:size]+"C"*(size-8)
    seq=major_allele
  
    seq={key:transform_from_fasta_to_index(key,ma) for key,ma in zip(seq,major_allele)}
    lik=generate_likelihood(trans_gen, ems_gen, initial, seq, char_map)
    return lik
    
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
    return band



def find_smallet_equivalence_class(params):
    size=2
    best_combi=params
    while size<=len(params):
        best_combi=get_best_of_size(best_combi, size)
        size*=2
    return best_combi
    
        


if __name__ == '__main__':
    #call_helper()
    print(call_forward())
    #lik=test_model_likelihood(8)
    #print(maximize_likelihood_exhaustive(lik,["pop1","pop2"]))
    print(find_smallet_equivalence_class(["a","b","a","a","a","b","a","a"]))
    
    