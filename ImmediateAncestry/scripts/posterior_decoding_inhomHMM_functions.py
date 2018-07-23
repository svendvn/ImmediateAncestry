import tmhmm
import numpy
import operator
import functools
import matrix_generation
import operator
import matplotlib.pyplot as plt
#from plot_chromosome import plot_chromosome
from math import log
from itertools import product
from argparse import ArgumentParser
from ggrandparents_model import generate_likelihood_from_data
from simulate_data import simulate, simulate_recombs, simulate_allele_frequencies
from recombination import read_recombination_file
from construct_seqs import get_seqs
from brute_force_maximization import maximize_likelihood_exhaustive
from likelihood_evaluations import evaluate_list, parse_configs
from shortcut_names import read_shortcuts
from setup_MCMC import mcmc_search
from id_dic import id_dic
from mock_likelihood import mock_likelihood
from matrix_generation import generate_emission_matrix
from matrix_generation import generate_transition_matrix
from ggrandparents_model import find_smallet_equivalence_class
from New_transition_emission_Matrix import initial_matrix
from matrix_generation import generate_transition_matrix


code_to_index={"A":0, "C":1, "G":2, "T":3}

# -----------------------------------------
# -----------------------------------------
# -----------------------------------------
def read_config_file (name, generations):
	filepath="../../../prepared_data/configs_Likelihood/"+name+"_"+str(generations)+".out"
	#filepath=name
	cols=2*2**generations
	config=[]
	file_object  = open(filepath, "r")
	for conf in file_object.read(cols):
		if conf != " ":
			config.append(conf)
	return config

def extrat_conf_index_from_emMat(ancestral_allele_dictionary,s,full_to_short,max_numberof_recombs,num_recombs):
    m=2**(s-1)
    M=m**2
    states=[]
    def generate(params):
        pops_mother=params[:m]
        pops_father=params[m:]
	 
        for i in range(m):   #the m ancestries of the mother
            for j in range(m): #the m ancestries of the father
                i1=i*m+j
                mother=full_to_short[pops_mother[i]]
                father=full_to_short[pops_father[j]]
                estado=mother+""+father
                states.append(estado)

        st=[]
        for i in (states):
            st.append(i)

        if max_numberof_recombs != -1:
            for i in range (num_recombs):
                for j in st:
                    states.append(j)	

        return states
    return generate




def best_path(Post_decod_matrix):
	keys=[]	#array of keys ( smllest configs)
	configs={}	#dic {'ee'=[[pos1, pos4],[probIn1,probIn4]]}
	probs=[]

	emptylist=[]

	for key1 in Post_decod_matrix:
		keys.append(key1)
		configs[key1]=[[],[]]

	for i in range(len(Post_decod_matrix[keys[0]])):
		prob1=0
		conf=''
		lastkey=''
		for key in keys:
			prob2=Post_decod_matrix[key][i]				

			if prob1 < prob2:
				prob1=prob2
				conf=key
			elif prob1 == prob2 :
				if key == lastkey:
					prob1=prob2
					conf=key		 
		superlista=configs[conf]
		lastkey=conf
		pos=superlista[0]
		pos.append(i)
		prob=superlista[1]
		prob.append(prob1)
		configs[conf]=[pos,prob]

	return configs


def calc_colors_PostDecod(best_config):
	color_dictionary={}
	color=0
	for key in best_config:
		pos=best_config[key][0]
		prob=best_config[key][1]
		print ("\nkey", key, "pos",len(pos))
		color=color+1
		color_dictionary[key]=color
	return color_dictionary



def posterior_decoding(alleles_list, recomb_map_list, seq_list, generations,max_numberof_recombs, num_recombs, params,full_to_short=id_dic(), short_to_full=id_dic()):

	#chrmos=['chr1','chrm2a','chrm2b','chrm3','chrm4','chrm5','chrm6','chrm7','chrm8','chrm9','chrm10','chrm11','chrm12','chrm13','chrm14','chrm15','chrm16','chrm17','chrm18','chrm19','chrm20','chrm21','chrm22']

	Post_dec_best_configs=[]
	for alleles, recomb_map, seq in zip(alleles_list, recomb_map_list, seq_list):

		#postFile="posteriordecoding_"+chrma+".txt"
		#f=open(postFile,"w")

		trans_gen=generate_transition_matrix(recomb_map, generations, False, max_numberof_recombs, num_recombs)

		ems_gen=matrix_generation.generate_emission_matrix(alleles, generations,max_numberof_recombs, num_recombs)
		params2=[short_to_full[p] for p in params]
		emission_generator=ems_gen(params2)

		states=extrat_conf_index_from_emMat(alleles, generations,full_to_short,max_numberof_recombs, num_recombs)
		initial=initial_matrix(generations,max_numberof_recombs,num_recombs)

		
		char_map={str(i):i for i in range(6)}
	
		seq2="".join(map(str,seq))

		print("Trans_Mat\n", trans_gen(0))
		print("Ems_Mat\n", emission_generator(0))
		print("Initial\n", initial)


		fv, constant= tmhmm.hmm.forward2(seq2, numpy.array(initial), trans_gen, emission_generator, char_map, None, None)

		bv = tmhmm.hmm.backward2(seq2, constant,numpy.array(initial), trans_gen, emission_generator, char_map, None, None)
		
		print("forward \n",fv)
		print("backward \n",bv)

		matrix=numpy.multiply(fv,bv)
		print ("matrix \n", matrix)

		Cj=numpy.sum(matrix, axis=1)
		print ("Cj\n",Cj)
	
		Cjreshape=Cj.reshape(len(Cj),1)
		print("Cjreshape\n", Cjreshape)

		Pmatrix=numpy.divide(matrix,Cjreshape)
		print("Pmatrix\n", Pmatrix)
	
		states_func=extrat_conf_index_from_emMat(alleles,generations,full_to_short,max_numberof_recombs,num_recombs)
		states=states_func(params2)
		correct_states=[]
		for i in states:
			correct_states.append(find_smallet_equivalence_class(i)) #Configs in small eq class
	
		print("correct_states\n",correct_states)

		#f.writelines(["%s\t" % item  for item in correct_states])
		
		#f.writelines(["%s\n" % line  for line in Pmatrix])
		#f.close() 
			#print("Pmatrix shape",Pmatrix.shape)
		#print("Pmatrix size", Pmatrix.size)
		Post_decod_dic={} #Dictionary with prob for each config
		for i in range(len(correct_states)):
			#print(correct_states[i])
	
			if correct_states[i] in Post_decod_dic:
			
				input_val=[Post_decod_dic[correct_states[i]],Pmatrix[:,i]]
				#print(type(input_val))
				#Post_decod_dic[correct_states[i]]=(sum(x) for x in zip(*input_val))
				Post_decod_dic[correct_states[i]]=list(map(sum, zip(*input_val)))
				#print(Post_decod_dic[correct_states[i]])
			else :
				Post_decod_dic[correct_states[i]]=Pmatrix[:,i]
	

		best_config = best_path(Post_decod_dic)

		color_dictionary=calc_colors_PostDecod(best_config)
			
		Post_dec_best_configs.append(dict(best_config))

		#if chrmos == 3:
		#	f.write("Chromosom: ")
		#	f.write(str(chrmos))
		#	f.write(str(best_config))
		#chrmos=chrmos+1

		
	return Post_dec_best_configs,color_dictionary













