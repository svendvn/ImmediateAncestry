from matrix_generation import generate_emission_matrix, generate_transition_matrix, initial_matrix
import tmhmm
import numpy

class id_dic(object):
    def __getitem__(self, key):
        return key
    
def plot_hidden_states(params, hidden_states):
    
    

def sim_hidden_states(alleles_list, recomb_map_list, params, seq_list, generations,max_numberof_recombs,num_recombs ,full_to_short=id_dic(), short_to_full=id_dic()):
    #f=open("filepos.txt","w") 
    #chrmos=1
    #chrmos=['chrm1','chrm2a','chrm2b','chrm3','chrm4','chrm5','chrm6','chrm7','chrm8','chrm9','chrm10','chrm11','chrm12','chrm13','chrm14','chrm15','chrm16','chrm17','chrm18','chrm19','chrm20','chrm21','chrm22']
    hidden_states=[]
    for alleles, recomb_map, seq in zip(alleles_list, recomb_map_list, seq_list):

        #postFile="simHD_"+chrma+".txt"
        #f=open(postFile,"w")

        trans_gen=generate_transition_matrix(recomb_map, generations, False, max_numberof_recombs, num_recombs)

        ems_gen=generate_emission_matrix(alleles, generations,max_numberof_recombs, num_recombs)

        params2=[short_to_full[p] for p in params]
        emission_generator=ems_gen(params2)

        #extrat_conf_index_from_emMat(alleles,generations,full_to_short,max_numberof_recombs,num_recombs)

        initial=initial_matrix(generations)
        
        char_map={str(i):i for i in range(6)}
    
        seq2="".join(map(str,seq))

    #    print("Trans_Mat\n", trans_gen(0))
    #    print("Ems_Mat\n", emission_generator(0))
    #    print("Initial\n", initial)

        fv, constant= tmhmm.hmm.forward2(seq2, numpy.array(initial), trans_gen, emission_generator, char_map, None, None)

        bv = tmhmm.hmm.backward2(seq2, constant,numpy.array(initial), trans_gen, emission_generator, char_map, None, None)

        matrix=numpy.multiply(fv,bv)
        Cj=numpy.sum(matrix, axis=1)
        Cjreshape=Cj.reshape(len(Cj),1)
        Pmatrix=numpy.divide(matrix,Cjreshape)
        
        xs=[]
        N,M=Pmatrix.shape        

        st=numpy.random.choice(list(range(M)),1, p=Pmatrix[0,:])

        xs.append(st)

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

        for j in range(1,len(seq)):
            
            em=emission_generator(j)
            tr=trans_gen(j-1)

        #    P=numpy.array(numpy.asmatrix(bv)[j,:] * em[:,seq[j]].T * tr[xs[j-1],:])
            P=numpy.array( bv[j,:] *numpy.squeeze(numpy.asarray(em[:,seq[j]].T))* numpy.squeeze(numpy.asarray(tr[xs[j-1],:])))    


            P=P/(numpy.sum(P))

            
            st=numpy.random.choice(list(range(M)),1, p=P)

            xs.append(st[0])
        hidden_states.append(xs)
    return hidden_states

