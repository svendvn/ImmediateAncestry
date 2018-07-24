from argparse import ArgumentParser
from numpy.random import choice
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

usage="""This program generates a likelihood given the inputs: sequence, allele frequencies, and genetic, recombinational distances.

It then searches the spaces of possible ancestries of the great grandparents of the input sequence.

The data pipeline is used for testing purposes because the program contains ways to simulate data according to the model. The default pipeline [6,7] will run the program as expected from the README.md.
"""

parser = ArgumentParser(usage=usage)

#input/output arguments

# # #
parser.add_argument("--seq_files", type=str, nargs="+", default=['../../../prepared_data/seqs/seq1.txt',
'../../../prepared_data/seqs/seq2a.txt','../../../prepared_data/seqs/seq2b.txt','../../../prepared_data/seqs/seq3.txt',
'../../../prepared_data/seqs/seq4.txt','../../../prepared_data/seqs/seq5.txt','../../../prepared_data/seqs/seq6.txt',
'../../../prepared_data/seqs/seq7.txt','../../../prepared_data/seqs/seq8.txt','../../../prepared_data/seqs/seq9.txt',
'../../../prepared_data/seqs/seq10.txt','../../../prepared_data/seqs/seq11.txt','../../../prepared_data/seqs/seq12.txt',
'../../../prepared_data/seqs/seq13.txt','../../../prepared_data/seqs/seq14.txt','../../../prepared_data/seqs/seq15.txt',
'../../../prepared_data/seqs/seq16.txt','../../../prepared_data/seqs/seq17.txt','../../../prepared_data/seqs/seq18.txt',
'../../../prepared_data/seqs/seq19.txt','../../../prepared_data/seqs/seq20.txt','../../../prepared_data/seqs/seq21.txt',
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


# # # Recombination Arguments 
parser.add_argument('--max_numberof_recombs', type=str, nargs='+', default=[-1], help="Argumnt to fix a maximun number of recombinations allowed per chromosome. Default is -1, which means that all the recombinations are allowed.")

parser.add_argument('--num_recombs', type=str, nargs='+', default=[2], help="Numer of recombinations allwed in case that the previous argument is enabled.")

# # # #

# # # Generation  and Index
parser.add_argument('--generations', type=int, default=2, help='The number of generations to go back. The number of ancestors to estimate will be 2 to the power of this number.')

parser.add_argument('--seq_indices', type=int, default=[0], nargs='+', help='the same as seq but here using indices instead of names.')
# # #



parser.add_argument('--outfiles', type=str, nargs='+', default=["immediate_ancestry_results.txt"], help="These are the files where the output will be written.")


#technical arguments
parser.add_argument("--seq", type=str, default=[], nargs="+", help="If seq_files contains many lines of corresponding to different ancestors, this list specifies which individual(s) to include in the analysis")



parser.add_argument('--ploidy_discrepancy', type=int, default=1, help='The seq_files are assumed to be diploid, but if they are haploid (and have no missing data) this can be set to 2, such that the program will combine the two specified haploid sequences into a diploid sequence.')

parser.add_argument("--pops_to_search", type=str, nargs="+",default=[], help="DEPRECATED. This is a list of populations to search for the best solution. Defaults to all available population.")

parser.add_argument('--auto_outfile_name', action='store_true', default=False, help="This will construct a filename from the names within the seq files and append the outfile-string and then '.txt'.")

parser.add_argument('--auto_outfile_id', type=str, default='0', help='id is an identifier for the output file.')

parser.add_argument('--truncate_af', type=float, default=0.01, help='This will truncate the allele frequencies such that the lowest value is truncate_af and the highest value is 1-truncate_af')

parser.add_argument('--mock_likelihood', action='store_true', default=False, help='For testing. this will replace the costly likelihood with a likelihood that is just minus the hamming distance.')

#simulation arguments
parser.add_argument('--sequences_pipeline', type=int, nargs='+', default=[6,7], help='This is the list of data processing steps to go through. 1 is from scratch, no inputs. 2 is deprecated. 3 is only allele frequencies. 4 is allele frequencies and ancestral sequences. 5 is allele frequencies, ancestral sequences and recombination. 6 is allelfreqs, recombs and sequences. 7 is the same as 6 but packed nicely together.')

parser.add_argument('--simulate_recombinations', default=False, action='store_true', help='If one doesnt provide a recombination map, the program will throw an error unless this is flag is turned on')

parser.add_argument('--sim_recomb_rate', type=float, default=0.1, help="If recombination rate is simulated, this is the average recombination probability.")
parser.add_argument('--recomb_rate_infinity', default=False, action='store_true', help='This will set the recombination rate to infinity, meaning that there is assumed no linkage.')
 
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

#annealing arguments

options = parser.parse_args()

if options.shortcut_names:

    full_to_short,short_to_full=read_shortcuts(options.shortcut_names)

else:
    class id_dic(object):
        def __getitem__(self, key):
            return key
    full_to_short, short_to_full=id_dic(), id_dic()
    
options.short_to_full=short_to_full
    
if options.true_pops:
    print(options.true_pops)

    options.true_pops=parse_configs(options.true_pops, options.generations, short_to_full=id_dic())
    print(options.true_pops)

if options.simulate_recombinations:
    pass # here we pass because this means that there is no information about the number of SNPs based on the recombination map recombs=get_recombinations(options.recomb_map, )
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
    
#    print(options.outfiles)
#    print(len(sequences))
            
    assert len(sequences)==len(options.outfiles)
    assert set(full_to_short.keys())==set(extra_info['pop_names'])

# -------------------------           

likelihoods=[generate_likelihood_from_data(allele_frequencies, 
                                           recombs, seq_system,  options.max_numberof_recombs,noptions.um_recombs ,
                                           options.generations, 
                                           short_to_full,
                                           rho_infinity=options.recomb_rate_infinity) for seq_system in sequences]
likelihood=likelihoods[0]

#print(likelihood.sequence)

#import sys
#sys.exit()
#print('likelihoods', likelihoods)

if options.mock_likelihood:
    likelihoods=[mock_likelihood for _ in likelihoods]
if options.type_of_analysis=='brute-force':
    popn=[full_to_short[n] for n in extra_info['pop_names']]
    ad=[]
    for likelihood in likelihoods:
        ad.append(maximize_likelihood_exhaustive(likelihood,  popn, options.generations))

if options.type_of_analysis=='mcmc_search':
    print ("options.type_of_analysis=='mcmc_search'")
    popn=[full_to_short[n] for n in extra_info['pop_names']]
    print ("popn", popn)
    ad=[]
    for likelihood in likelihoods:
        ad.append(mcmc_search(likelihood, recombs, sequences, popn, options.generations, short_to_full, N=options.mcmc_reps))
        print ("ad", ad)

res=''

for bd,outfile in zip(ad, options.outfiles):
    # bd :  [('ettv', -13093.757802897619), ('tvtv', -13094.839761384726),....
    # outfile : 
    res=''
    for params,likel in bd:
        print("params = ", params)
        print("likel = ", likel)
        full_names=list(map(str,params))
        print ("full_names = ", full_names)
        short_names=[full_to_short[name] if len(name)>1 else name for name in full_names]
        print ("short_names", short_names)
        res+=" ".join(short_names)+" "+str(likel)+"\n"
    print("res = ", res)
    with open(outfile   , "w") as f:
        f.write(res)

