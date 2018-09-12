from matrix_generation import generate_emission_matrix, generate_transition_matrix, initial_matrix
import tmhmm
import numpy
from collections import Counter
import subprocess
from six import string_types
import os


class id_dic(object):
    def __getitem__(self, key):
        return key
    
def plot_hidden_states(params, hidden_states, chrom_names, individual_name, shortcut_name_file, plot_filename=None):
    if not isinstance(params, string_types):
        params=''.join(params)
    if plot_filename is None:
        plot_filename=individual_name+'_'+params+'.pdf'
    temp_data_frame_name='.'.join(plot_filename.split('.')[:-1])+'.tmpcsv'
    with open(temp_data_frame_name, 'w') as f:
        for hidden_sequence, chrom_name in zip(hidden_states, chrom_names):
            for hidden_state in hidden_sequence:
                f.write(str(hidden_state)+','+str(chrom_name)+'\n')
    project_dir = os.path.dirname(__file__)
    command=['Rscript',os.path.join(project_dir, 'plot_hidden_states.R'), 
             temp_data_frame_name, 
             params, 
             shortcut_name_file, 
             individual_name, 
             plot_filename]
    subprocess.call(command)
    
def get_states(num, m):
    m=len(params)//2
    mother_i=num%m
    father_i=num//m
    return mother_i,  father_i

def refine_estimate(params, hidden_states, params):
    count_dic=Counter(params)
    m=len(params)//2
    mother_counts=[0]*m
    father_counts=[0]*m
    for double_states_index,  count_val in count_dic.items():
        mother_i,  father_i=get_states(double_states_index,  m)
        mother_counts[mother_i]+=count_val
        father_counts[father_i]+=father_val
    first_half=[count_dic[i] for i in range(len(params)//2)]
    second_half=[count_dic[i] for in range(len(params)//2,  len(params))]
    
    for hidden_sequence in hidden_states:
        for hidden_state in hidden_sequence:
            get_states()
            
    
    
    
    

def sim_hidden_states(likelihoods, params):
    #f=open("filepos.txt","w") 
    #chrmos=1
    #chrmos=['chrm1','chrm2a','chrm2b','chrm3','chrm4','chrm5','chrm6','chrm7','chrm8','chrm9','chrm10','chrm11','chrm12','chrm13','chrm14','chrm15','chrm16','chrm17','chrm18','chrm19','chrm20','chrm21','chrm22']
    hidden_states=[]
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
        
        xs=[]
        N,M=Pmatrix.shape        

        st=numpy.random.choice(list(range(M)),1, p=Pmatrix[0,:])

        xs.append(st[0])

        #states_func=extrat_conf_index_from_emMat(alleles, generations, full_to_short, max_numberof_recombs, num_recombs)
        #states=states_func(params2)

        #correct_states=[]
        #for i in states:
        #    correct_states.append(find_smallet_equivalence_class(i)) 
        #print ("correctstates", correct_states)
        #print ("Primer st", st)
        #print ("Primer correct state", correct_states[st])
        #    colores[i]=""
        
        #couples={}
        #couples[correct_states[st]]=[0]

        #f.write(str(correct_states[st]))
        #f.write(" ")

        for j in range(1,len(seq2)):
            
            log_em=emission_generator(j)
            tr=trans_gen(j-1)

        #    P=numpy.array(numpy.asmatrix(bv)[j,:] * em[:,seq[j]].T * tr[xs[j-1],:])
            b=bv[j,:]
            e=log_em[:,int(seq2[j])]
            t=numpy.log(tr[xs[j-1],:])
            #print(b,e,t)
            P=b+e+t 
            #print(P)
            
            P=P-numpy.max(P)
            #print(P)
            P=numpy.exp(P)
            #print(P)
            P=P/numpy.sum(P)
            #print(P)

            
            st=numpy.random.choice(list(range(M)),1, p=P)

            xs.append(st[0])
        hidden_states.append(xs)
    return hidden_states

