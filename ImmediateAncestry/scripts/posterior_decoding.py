import numpy
from copy import deepcopy
import tmhmm

def get_states(num, m):
    mother_i=num//m
    father_i=num%m
    return mother_i,  father_i

def refine_estimates_from_posterior_decoding(params, posterior_decodings):
    new_params=list(deepcopy(params))
    m=len(params)//2
    count_vector=[0]*(m**2)
    for posterior_decoding_matrix in posterior_decodings:
        new_v=numpy.sum(posterior_decoding_matrix,  axis=0)
        count_vector=[c+v for c, v in zip(count_vector,  new_v)]
    print('joint counts',count_vector)
    mother_counts=[0]*m
    father_counts=[0]*m
    for double_states_index,  count_val in enumerate(count_vector):
        mother_i,  father_i=get_states(double_states_index,  m)
        mother_counts[mother_i]+=count_val
        father_counts[father_i]+=count_val
    print('mother counts', mother_counts)
    print('father counts', father_counts)
    mimin, mmin=numpy.argmin(mother_counts),  numpy.min(mother_counts)
    mimax, mmax=numpy.argmax(mother_counts),  numpy.max(mother_counts)
    fimin, fmin=numpy.argmin(father_counts),  numpy.min(father_counts)
    fimax, fmax=numpy.argmax(father_counts),  numpy.max(father_counts)
    if mmin*3<mmax:
        print('exchanging a',  new_params[mimin],  'for a',  new_params[mimax])
        new_params[mimin]=new_params[mimax]
    if fmin*3<fmax:
        new_params[m+fimin]=new_params[m+fimax]
        print('exchanging a',  new_params[fimin],  'for a',  new_params[fimax])
    return ''.join(new_params)
    

def refine_estimate(likelihoods,  params):
    posterior_decodings=posterior_decoding(likelihoods,  params)
    return refine_estimates_from_posterior_decoding(params,  posterior_decodings)

def posterior_decoding(likelihoods, params):

	#chrmos=['chr1','chrm2a','chrm2b','chrm3','chrm4','chrm5','chrm6','chrm7','chrm8','chrm9','chrm10','chrm11','chrm12','chrm13','chrm14','chrm15','chrm16','chrm17','chrm18','chrm19','chrm20','chrm21','chrm22']

    Pmatrices=[]
    for likelihood in likelihoods:
        
        #postFile="simHD_"+chrma+".txt"
        #f=open(postFile,"w")

        trans_gen=likelihood.transition_generator

        ems_gen=likelihood.emission_generator_function

        params2=[likelihood.short_to_full[p] for p in params]
        emission_generator=ems_gen(params2, log=True)

        #extrat_conf_index_from_emMat(alleles,generations,full_to_short,max_numberof_recombs,num_recombs)

        initial=likelihood.initial_probabilities
        
        char_map=likelihood.char_map
    
        seq2=likelihood.sequence

    #    print("Trans_Mat\n", trans_gen(0))
    #    print("Ems_Mat\n", emission_generator(0))
    #    print("Initial\n", initial)

        fv, constant= tmhmm.hmm.forward2log(seq2, numpy.array(initial), trans_gen, emission_generator, char_map, None, None)
        #print('fv',fv)
        bv = tmhmm.hmm.backward2log(seq2, constant, numpy.array(initial), trans_gen, emission_generator, char_map, None, None)
        #print('bv',bv)
        matrix=fv+bv
        maxes=numpy.max(matrix, axis=1)
        matrix-=maxes[:,numpy.newaxis]
        matrix=numpy.exp(matrix)
        rowsums=numpy.sum(matrix, axis=1)
        Pmatrix=matrix/rowsums[:,numpy.newaxis]
        Pmatrices.append(Pmatrix)
    
    return Pmatrices
    
if __name__=='__main__':
    posterd=numpy.array([[0.4, 0.1, 0.4, 0.1],[0.4, 0.1, 0.4, 0.1], [0.2, 0.4, 0.02, 0.2, 0.4] ])
    print(posterd)
    print(refine_estimate('tsev',  [posterd]))
    
