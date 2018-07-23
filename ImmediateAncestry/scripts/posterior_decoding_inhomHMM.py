import tmhmm
import numpy
import operator
import functools
import matrix_generation
import operator
import matplotlib.pyplot as plt
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

from posterior_decoding_inhomHMM_functions import posterior_decoding
from posterior_decoding_inhomHMM_functions import best_path
from posterior_decoding_inhomHMM_functions import read_config_file
from posterior_decoding_inhomHMM_functions import extrat_conf_index_from_emMat


code_to_index={"A":0, "C":1, "G":2, "T":3}

# -----------------------------------------
# -----------------------------------------
# -----------------------------------------

usage="""This program generates a posterior decoding matrix and a kariotype plot showing the most probably configuration of ancestors at each position. Mandatory inputs are: sequence, allele frequencies, and genetic, recombinational distances."""
parser = ArgumentParser(usage=usage)

# # # #
parser.add_argument("--seq_files", type=str, nargs="+", 
default=['../../../prepared_data/seqs/seq1.txt','../../../prepared_data/seqs/seq2a.txt',
'../../../prepared_data/seqs/seq2b.txt','../../../prepared_data/seqs/seq3.txt',
'../../../prepared_data/seqs/seq4.txt','../../../prepared_data/seqs/seq5.txt',
'../../../prepared_data/seqs/seq6.txt','../../../prepared_data/seqs/seq7.txt',
'../../../prepared_data/seqs/seq8.txt','../../../prepared_data/seqs/seq9.txt',
'../../../prepared_data/seqs/seq10.txt','../../../prepared_data/seqs/seq11.txt',
'../../../prepared_data/seqs/seq12.txt','../../../prepared_data/seqs/seq13.txt',
'../../../prepared_data/seqs/seq14.txt','../../../prepared_data/seqs/seq15.txt',
'../../../prepared_data/seqs/seq16.txt','../../../prepared_data/seqs/seq17.txt',
'../../../prepared_data/seqs/seq18.txt','../../../prepared_data/seqs/seq19.txt',
'../../../prepared_data/seqs/seq20.txt','../../../prepared_data/seqs/seq21.txt',
'../../../prepared_data/seqs/seq22.txt'], help="This is a list of files. Each file correspond to a different chromosome. The files should contain space separated lines of allele types. 0,1,2 means that there are 0,1 or 2 copies of the allele. 3 is missing data.")

parser.add_argument("--recomb_map", type=str, default=['../../../prepared_data/rhos/rho_chr1.txt',
'../../../prepared_data/rhos/rho_chr2a.txt','../../../prepared_data/rhos/rho_chr2b.txt',
'../../../prepared_data/rhos/rho_chr3.txt','../../../prepared_data/rhos/rho_chr4.txt',
'../../../prepared_data/rhos/rho_chr5.txt','../../../prepared_data/rhos/rho_chr6.txt',
'../../../prepared_data/rhos/rho_chr7.txt','../../../prepared_data/rhos/rho_chr8.txt',
'../../../prepared_data/rhos/rho_chr9.txt','../../../prepared_data/rhos/rho_chr10.txt',
'../../../prepared_data/rhos/rho_chr11.txt','../../../prepared_data/rhos/rho_chr12.txt',
'../../../prepared_data/rhos/rho_chr13.txt','../../../prepared_data/rhos/rho_chr14.txt',
'../../../prepared_data/rhos/rho_chr15.txt','../../../prepared_data/rhos/rho_chr16.txt',
'../../../prepared_data/rhos/rho_chr17.txt','../../../prepared_data/rhos/rho_chr18.txt',
'../../../prepared_data/rhos/rho_chr19.txt','../../../prepared_data/rhos/rho_chr20.txt',
'../../../prepared_data/rhos/rho_chr21.txt','../../../prepared_data/rhos/rho_chr22.txt'], nargs="+", help="This is a list of files containing the genetic distances between SNPs. Each file correspond to one chromosome and contains a single space-separated line of N-1 genetic distances, where N is the number of SNPs. The genetic distance is measured in probability of recombination in one generation between the two SNPs.")

parser.add_argument("--allele_frequencies", type=str, default=['../../../prepared_data/freqs/chr1_freqs.txt',
'../../../prepared_data/freqs/chr2a_freqs.txt','../../../prepared_data/freqs/chr2b_freqs.txt',
'../../../prepared_data/freqs/chr3_freqs.txt','../../../prepared_data/freqs/chr4_freqs.txt',
'../../../prepared_data/freqs/chr5_freqs.txt','../../../prepared_data/freqs/chr6_freqs.txt',
'../../../prepared_data/freqs/chr7_freqs.txt','../../../prepared_data/freqs/chr8_freqs.txt',
'../../../prepared_data/freqs/chr9_freqs.txt','../../../prepared_data/freqs/chr10_freqs.txt',
'../../../prepared_data/freqs/chr11_freqs.txt','../../../prepared_data/freqs/chr12_freqs.txt',
'../../../prepared_data/freqs/chr13_freqs.txt','../../../prepared_data/freqs/chr14_freqs.txt',
'../../../prepared_data/freqs/chr15_freqs.txt','../../../prepared_data/freqs/chr16_freqs.txt',
'../../../prepared_data/freqs/chr17_freqs.txt','../../../prepared_data/freqs/chr18_freqs.txt',
'../../../prepared_data/freqs/chr19_freqs.txt','../../../prepared_data/freqs/chr20_freqs.txt',
'../../../prepared_data/freqs/chr21_freqs.txt','../../../prepared_data/freqs/chr22_freqs.txt'], nargs="+", help="This is a list of files containing the allele frequencies in the known populations. Each file correspond to a chromosome and they contain space_separated lines of the allele frequencies of length N genes.")
# # # #

# # # Recombination Arguments 
parser.add_argument('--max_numberof_recombs', type=str, nargs='+', default=1, help="Argument to set a maximun number of recombinations allowed per chromosome. Default is -1, which means that infinite recombinations are allowed.")

parser.add_argument('--num_recombs', type=str, nargs='+', default=1, help="Numer of recombinations allwed in case that the previous argument is enabled.")




# # # #

# # # Generations and Index
parser.add_argument('--generations', type=int, default=2, help='The number of generations to go back. The number of ancestors to estimate will be 2 to the power of this number.')

parser.add_argument('--seq_indices', type=int, default=[218], nargs='+', help='the same as seq but here using indices instead of names.')
# # #



# # # # Test
#parser.add_argument("--seq_files", type=str, nargs="+", default=['../../../prepared_data/Test_Files/seq.txt'], help="This is a list of files. Each file correspond to a different chromosome. The files should contain space separated lines of allele types. 0,1,2 means that there are 0,1 or 2 copies of the allele. 3 is missing data.")

#parser.add_argument("--recomb_map", type=str, default=['../../../prepared_data/Test_Files/rhos.txt'], nargs="+", help="This is a list of files containing the genetic distances between SNPs. Each file correspond to one chromosome and contains a single space-separated line of N-1 genetic distances, where N is the number of SNPs. The genetic distance is measured in probability of recombination in one generation between the two SNPs.")

#parser.add_argument("--allele_frequencies", type=str, default=['../../../prepared_data/Test_Files/freqs.txt'], nargs="+", help="This is a list of files containing the allele frequencies in the known populations. Each file correspond to a chromosome and they contain space_separated lines of the allele frequencies of length N genes.")
# # # #



parser.add_argument('--outfiles', type=str, nargs='+', default=["immediate_ancestry_results.txt"], help="These are the files where the output will be written.")
parser.add_argument("--seq", type=str, default=[], nargs="+", help="If seq_files contains many lines of corresponding to different ancestors, this list specifies which individual(s) to include in the analysis")
parser.add_argument('--ploidy_discrepancy', type=int, default=1, help='The seq_files are assumed to be diploid, but if they are haploid (and have no missing data) this can be set to 2, such that the program will combine the two specified haploid sequences into a diploid sequence.')
parser.add_argument("--pops_to_search", type=str, nargs="+",default=[], help="DEPRECATED. This is a list of populations to search for the best solution. Defaults to all available population.")
parser.add_argument('--auto_outfile_name', action='store_true', default=False, help="This will construct a filename from the names within the seq files and append the outfile-string and then '.txt'.")
parser.add_argument('--auto_outfile_id', type=str, default='0', help='id is an identifier for the output file.')
parser.add_argument('--truncate_af', type=float, default=0.01, help='This will truncate the allele frequencies such that the lowest value is truncate_af and the highest value is 1-truncate_af')
parser.add_argument('--mock_likelihood', action='store_true', default=False, help='For testing. this will replace the costly likelihood with a likelihood that is just minus the hamming distance.')
parser.add_argument('--sequences_pipeline', type=int, nargs='+', default=[6,7], help='This is the list of data processing steps to go through. 1 is from scratch, no inputs. 2 is deprecated. 3 is only allele frequencies. 4 is allele frequencies and ancestral sequences. 5 is allele frequencies, ancestral sequences and recombination. 6 is allelfreqs, recombs and sequences. 7 is the same as 6 but packed nicely together.')
parser.add_argument('--simulate_recombinations', default=False, action='store_true', help='If one doesnt provide a recombination map, the program will throw an error unless this is flag is turned on')
parser.add_argument('--recomb_rate', type=float, default=0.1, help="If recombination rate is simulated, this is the average recombination probability.")
parser.add_argument('--skewness', type=float, default=4, help="The skewness of the simulated recombination rates. If set to 0, there will be a constant recombination rate. Values below -1 are not meaningful")
parser.add_argument('--ancestors', type=str, default=[], nargs='+', help="This is the names of the ancestors in the ancestor files which are chosen as ancestors if the sequences are simulated from non-simulated ancestors.")
parser.add_argument("--ancestor_indices", type=int, default=[], nargs="+", help="This list of numbers are the indices of the ancestors in the ancestor files which are chosen as ancestors.")
parser.add_argument("--ancestor_files", type=str, default=[], nargs="+", help="This is a list of files containing the ancestor alleles. The setup is the same as the seq_files. It should only be applied if sequences should be simulated.")
parser.add_argument("--SNPs", type=int, default=0, help="this is the number of SNPs to simulate per segment. If any input files are specified this will be ignored(because then the number of SNPs will match the number of SNPs in the input files).")
parser.add_argument("--mcmc_reps", type=int, default=300, help="this is the number of runs of the MCMC (if used).")
parser.add_argument("--reps", type=int, default=1, help="This is the number of independent segments of SNPs to draw if simulated")
parser.add_argument("--true_pops", type=str, nargs="+", default=[], help="If simulations take place, this will be the true ancestors. It has to be of length 2**generations*ancestor_multiplier")
parser.add_argument('--no_pops', type=int, default=4, help='the number of populations.')
parser.add_argument('--ancestor_multiplier', type=int, default=1, help='the number of sets of ancestors are simulated. If more than one, the analysis will also be run more than once.')
parser.add_argument('--sequence_multiplier', type=int, default=1, help='the number of sequences that should be simulated from each set of ancestors.')
parser.add_argument('--population_names', type=str, nargs='+', default=[], help='When allele frequencies are simulated, the population names can not be read from anywhere, so to avoid default values (a,b,c,d, etc.) it can be specified here.')
parser.add_argument('--type_of_analysis', type=str, choices=['brute-force', 
                                                             'mcmc_search', 
                                                             'evaluate_likelihoods'], 
                    default='mcmc_search',
		#default='brute-force',
                    help='chooses the type of analysis to run on the data. \
                          brute-force searches all possible combinations of gparents (based on the populations specified in allele_frequencies(that may be simulated)). It is strongly adviced not to use this setting for generations>=4. \
                          evaluate_likelihoods evaluates the likelihoods specified in the list --configs_to_test \
                          mcmc_search is a false MCMC, that finds the maximum likelihood.')
parser.add_argument('--configs_to_test', type=str, nargs='+', default=['trivial_ellioti2.txt'], help='If type_of_analysis is evaluate_likelihoods this is a file(which has to contain a dot) of all the likelihoods to evaluate. It has to be space separated and a configuration on each line.')
parser.add_argument('--shortcut_names', type=str, default='shortcut_names2.txt', help='file that short cuts long names for easier readability. It is of the form [full_name short_name\n,...]')



options = parser.parse_args()
if options.shortcut_names:
    full_to_short,short_to_full=read_shortcuts(options.shortcut_names)
else:
    class id_dic(object):
        def __getitem__(self, key):
            return key
    full_to_short, short_to_full=id_dic(), id_dic()
options.short_to_full=short_to_full
options.full_to_short=full_to_short
if options.true_pops:
    options.true_pops=parse_configs(options.true_pops, options.generations, short_to_full=id_dic())
if options.simulate_recombinations:
    pass # here we pass because this means that there is no information about the number of SNPs based on the 		recombination map recombs=get_recombinations(options.recomb_map, )
elif options.recomb_map:
    recombs, setups=read_recombination_file(options.recomb_map)
    extra_info={'setups':setups}
    extra_info={'recombs':recombs}

else:
    assert options.recomb_map, 'Recombination map not specified' 
outp = get_seqs(options,  extra_info)
sequences, allele_frequencies, recombination_map= outp
if options.auto_outfile_name:
    outfiles=[]
    for true_p in options.true_pops:
        for j in range(options.sequence_multiplier):
            filename=''.join(true_p)+'_'+str(j)+'_'+str(options.auto_outfile_id)
            outfiles.append(filename)
    options.outfiles=outfiles 
    assert len(sequences)==len(options.outfiles)
    assert set(full_to_short.keys())==set(extra_info['pop_names'])


# -------------------------------

def plot_postDecodLine (best_configs,color_dic):
	axisX=0
	num_colors=len(color_dic)
	alty=0


	for configs in best_configs:
		axisX=0
		alty=alty+1

		for key in configs:
			pos=configs[key][0]
			prob=configs[key][1]
			axisX=axisX + len(pos)
			
			plt.plot(pos,prob, "-",label=key)

		plt.axis([0,axisX , 0.4, 1.2])
		plt.legend(loc='best')
	#plt.plot(configs)
	#plt.ylabel('Probability')
	#plt.xlabel(paths)
		plt.show()

def plot_Chrom_simple (Post_dec_best_configs,color_dic):
	#print("POSTTTTTTTTTTTTTTTTTTTTTTT",len(Post_dec_best_configs))
	final_axisX=0
	alty=0
#	possible_colors = ['none','green','red','turquoise','orange','royalblue','blue','black','fuchsia',
#'yellow','salmon','sienna','chartreuse','lightseagreen','silver','darkcyan','blueviolet', 'purple','pink','tan', 'olivedrab','tomato','moccasin','gold','mediumspringgreen','brown1','cadetblue4','cadmiumorange','chartreuse1', 'cobalt','cobaltgreen','crimson','darkorchid4','deeppink4','gray59','hotpink','lightblue1','lightslateblue']

	possible_colors = {'ee':'green','et':'red','es':'turquoise','ev':'orange','tt':'royalblue', 'ts':'blue','tv':'black','ss':'fuchsia',
'sv':'yellow','vv':'pink','vs':'sienna','vt':'chartreuse','ve':'lightseagreen','st':'silver','se':'darkcyan', 'te':'blueviolet'}
	i=0
	for configs in Post_dec_best_configs:
		axisX=0
		alty=alty+1

		for key in configs:
			pos=configs[key][0]
	
			prob=configs[key][1]
			axisX=axisX + len(pos)
			y=[alty]*len(pos)
		#	plt.plot(pos,y, "|",label=key, markersize=5,color=possible_colors[color_dic[key]])
			plt.plot(pos,y, "|",label=key, markersize=5,color=possible_colors[key])
			if i<len(color_dic):
				plt.legend(loc='best', markerscale=0.5, scatterpoints=1, fontsize=10 ,numpoints=50)
				i=i+1

		if final_axisX< axisX:
			final_axisX= axisX
	#plt.legend(loc='best', markerscale=0.5, scatterpoints=1, fontsize=10 ,numpoints=10)	
	plt.axis([-10,final_axisX+100 , 0, alty+0.5])
	
	#plt.plot(configs)
	#plt.ylabel('Probability')
	#plt.xlabel(paths)   
#	plt.show()
	plt.savefig('PostDec_i218_Kidongo_2gen_rhoNoInf_1recMax.png')

# -------------------------------  int main  ------------------------------- #


#print("allele_frequencies\n",allele_frequencies)
#print("recombs\n",recombs)
#print("sequences\n",sequences)
#print("generations",options.generations)


#name=['i0_Brigitte','i1_11046','i2_Dandy','i3_11371','i4_11571','i5_11664','i6_Belle','i7_Shirley','i8_Micheline','i9_Paula','i10_11796','i11_Roosje','i12_11837','i13_11932','i14_Bonobo','i15_Julie','i16_Arnold','i17_Agathe','i18_12103','i19_Jannis','i20_12264_13649','i21_Chispi','i22_Gypso','i23_Giambo','i24_Paulinchen','i25_Fleur','i26_Peter','i27_12478','i28_Amadine','i29_Tarzan','i30_Patty','i31_Dirch','i32_Qafzeh','i33_12653','i34_12660','i35_12664','i36_Hanna','i37_Jeremy','i38_12695','i39_Arthur','i40_Cookie','i41_12771','i42_12788','i43_Josef','i44_12876','i45_Kindia','i46_Tommy','i47_12944','i48_12949','i49_Camilla','i50_13051','i51_Zieglein','i52_Knerten','i53_13123','i54_13125','i55_Negine','i56_Julius_jr','i57_13217','i58_13283','i59_Cathy','i60_13460','i61_13477','i62_Lucas','i63_Wingu','i64_Bingo','i65_Buddy','i66_Fiffy','i67_Peggy','i68_Wilma','i69_Achille','i70_Babsi','i71_Patrick','i72_Prudence','i73_13692','i74_Fritz','i75_Frederike','i76_Jonas','i77_13748','i78_Cindy','i79_Julchen','i80_Penny_Dicke','i81_13790','i82_13804','i83_Alisa','i84_Siggi','i85_13810','i86_Ronja','i87_13824','i88_Regina','i89_Amelie','i90_Loulou','i91_13942','i92_Babsy','i93_Bobo','i94_Edward','i95_Ezo','i96_Fanny','i97_Guille','i98_Ivan','i99_Jojo','i100_Judi','i101_Kyo','i102_Lola','i103_Louise','i104_Marlin','i105_Frits','i106_Louise','i107_Zwartje','i108_Sammy','i109_11664','i110_Belle','i111_Paula','i112_11815','i113_Benjie','i114_Bonobo','i116_Szymi','i117_12349','i118_Scholzi','i119_12396','i120_Fleur','i121_Mumin','i122_Amadine','i123_Tarzan','i124_Bingo','i125_Patty','i126_Willy','i127_Gipsy','i128_12591','i129_12604','i130_Anna_DE390','i131_Anna_S651','i132_Jeremy','i133_12695','i134_Julie','i135_Tommy','i136_Chura','i137_12962','i138_Lindi','i139_Camilla','i140_Karibuna','i141_13122','i142_13306','i143_Cathy','i144_13477','i145_13496','i146_13502','i147_Panya','i148_13709','i149_Johnny','i150_Konda','i151_Frederike','i152_13790','i153_Giuditta','i154_Churrero','i155_Marina','i156_Clara','i157_Silvia','i158_Chigo','i159_14088','i160_Marco','i161_14093','i162_Giudy','i163_Bongo','i164_Klaus','i165_Ali_Kaka','i166_Bashu','i167_Big_Jane','i168_Bili','i169_Cleo_DE068','i170_Cleo_S656','i171_Coco_D066','i172_Coco_S654','i173_Cumbo','i174_Diana_D059','i175_Diana_S655','i176_Doris','i177_Hans','i178_Jac','i179_Kambo','i180_Karla','i181_Kina','i182_Max','i183_Mgbadolite','i184_Mwisho','i185_Naika','i186_Natasha','i187_Noel','i188_Pipa','i189_Roy','i190_Safari','i191_Seki','i192_Socrates','i193_Sultana','i194_Tika','i195_Trixie_DE058','i196_Trixie_S657','i197_Yoyo','i198_Maki','i199_Manu','i200_Mgbadolite','i201_Diana','i202_Ikuru','i203_Trixie','i204_Akwaya_Jean','i205_Banyo','i206_Basho','i207_Damian','i208_Julie','i209_Kopongo','i210_Koto','i211_Paquita','i212_Taweh','i213_Tobi','i214_Vincent','i215_Andromeda','i216_Harriet','i217_Bwambale','i218_Kidongo','i219_Nakuu','i220_Padda','i221_Cindy','i222_Frederike','i223_Washu','i224_Athanga','i225_Coco','i226_Mgbadolite','i227_Tongo','i228_Cleo','i229_Bihati','i230_Maya','i231_Cindy','i232_Mirinda','i233_Alfred','i234_Ula','i235_Lara','i236_Luky','i237_Gamin','i238_Brigitte','i239_Vaillant','i240_Doris','i241_Julie','i242_Clara','i243_Noemie','i244_Yogui','i245_Tibe','i246_Blanquita','i247_Negrita','i248_Marlin','i249_Bosco','i250_Donald','i251_Jimmie','i252_Berta','i253_Annie','i254_Mike','i256_Linda','i259_Cindy','i261_Koby','i262_Peggy','i263_Salleh','i264_Toti','i265_14069','i266_Yaki']


#name="i250_Donald"
#name="i0_Brigitte"
#name="i242_Clara"
name="i218_Kidongo"

params=read_config_file(name, options.generations)

Post_best,color_dic=posterior_decoding(allele_frequencies, recombs, sequences[0],options.generations, options.max_numberof_recombs, options.num_recombs,params,full_to_short, short_to_full)

## -- Choose one of this repreentations:

#plot_postDecodLine(Post_best,color_dic)	# Posterior decoding probability fr each SNP
plot_Chrom_simple(Post_best,color_dic)	# Simple whoole chromosome with the configuration of each SNP





# ~~~~~~~~~~~~~~~~~~~~

#	for i in range ( #PAra ir cogiendo cada invividuo uno a uno
#	parser.add_argument('--seq_indices', type=int, default=[0], nargs='+', help='the same as seq but here using indices instead of names.')

#chromosomes =['1','2a','2b','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']







