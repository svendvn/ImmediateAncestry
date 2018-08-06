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
from print_structure import print_recombination_structure, print_sequence_structure

usage="""This program generates a likelihood given the inputs: sequence, allele frequencies, and genetic, recombinational distances.

It then searches the spaces of possible ancestries of the great grandparents of the input sequence.

The data pipeline is used for testing purposes because the program contains ways to simulate data according to the model. The default pipeline [6,7] will run the program as expected from the README.md.
"""

parser = ArgumentParser(usage=usage)

#input/output arguments
parser.add_argument("--seq_files", type=str, nargs="+", default=['/home/svendvn/Dropbox/Bioinformatik/Zoomodel/data/rawData2/hybrids/seq21.txt',
                                                                 '/home/svendvn/Dropbox/Bioinformatik/Zoomodel/data/rawData2/hybrids/seq22.txt'], help="This is a list of files. Each file correspond to a different chromosome. The files should contain space separated lines of allele types. 0,1,2 means that there are 0,1 or 2 copies of the allele. 3 is missing data.")
parser.add_argument("--recomb_map", type=str, default=['/home/svendvn/Dropbox/Bioinformatik/Zoomodel/data/rawData2/hybrids/rho_chr21.txt',
                                                       '/home/svendvn/Dropbox/Bioinformatik/Zoomodel/data/rawData2/hybrids/rho_chr22.txt'], nargs="+", help="This is a list of files containing the genetic distances between SNPs. Each file correspond to one chromosome and contains a single space-separated line of N-1 genetic distances, where N is the number of SNPs. The genetic distance is measured in probability of recombination in one generation between the two SNPs.")
parser.add_argument('--outfiles', type=str, nargs='+', default=["immediate_ancestry_results.txt"], help="These are the files where the output will be written.")
parser.add_argument("--allele_frequencies", type=str, default=['/home/svendvn/Dropbox/Bioinformatik/Zoomodel/data/rawData2/hybrids/chr21_freqs.txt',
                                                               '/home/svendvn/Dropbox/Bioinformatik/Zoomodel/data/rawData2/hybrids/chr22_freqs.txt'], nargs="+", help="This is a list of files containing the allele frequencies in the known populations. Each file correspond to a chromosome and they contain space_separated lines of the allele frequencies of length N genes.")

#model arguments
parser.add_argument('--generations', type=int, default=2, help='The number of generations to go back. The number of ancestors to estimate will be 2 to the power of this number.')

#technical arguments
parser.add_argument("--seq", type=str, default=[], nargs="+", help="If seq_files contains many lines of corresponding to different ancestors, this list specifies which individual(s) to include in the analysis")
parser.add_argument('--seq_indices', type=int, default=[0], nargs='+', help='the same as seq but here using indices instead of names.')
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
                    help='chooses the type of analysis to run on the data. \
                          brute-force searches all possible combinations of gparents (based on the populations specified in allele_frequencies(that may be simulated)). It is strongly adviced not to use this setting for generations>=4. \
                          evaluate_likelihoods evaluates the likelihoods specified in the list --configs_to_test \
                          mcmc_search is a false MCMC, that finds the maximum likelihood.')
parser.add_argument('--configs_to_test', type=str, nargs='+', default=['trivial_ellioti2.txt'], help='If type_of_analysis is evaluate_likelihoods this is a file(which has to contain a dot) of all the likelihoods to evaluate. It has to be space separated and a configuration on each line.')
parser.add_argument('--shortcut_names', type=str, default='shortcut_names2.txt', help='file that short cuts long names for easier readability. It is of the form [full_name short_name\n,...]')
parser.add_argument('--thin_coef', type=int, default=1, help='The thinning coefficient. If it is n, for each chromosome n new chromosomes will be made. If an original sequence is 1,2,3,4,5,.., a new sequences will be j,n+j,2n+j,.., for j=1,...,n. One has to put 7 8 in the sequences pipeline')
parser.add_argument('--')

#annealing arguments

options = parser.parse_args()

assert options.thin_coef==1 or 8 in options.sequences_pipeline, 'If thin_coef is not 1, 8 should be added to the covariance pipeline'

if options.shortcut_names:
    full_to_short,short_to_full=read_shortcuts(options.shortcut_names)
else:
    class id_dic(object):
        def __getitem__(self, key):
            return key
    full_to_short, short_to_full=id_dic(), id_dic()
    
options.short_to_full=short_to_full
    
if options.true_pops:
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

print_recombination_structure(recombination_map)
print_sequence_structure(sequences)

if options.auto_outfile_name:
    outfiles=[]
    for true_p in options.true_pops:
        for j in range(options.sequence_multiplier):
            filename=''.join(true_p)+'_'+str(j)+'_'+str(options.auto_outfile_id)
            outfiles.append(filename)
    options.outfiles=outfiles
    
    print(options.outfiles)
    print(len(sequences))
            
    assert len(sequences)==len(options.outfiles)
    
    assert set(full_to_short.keys())==set(extra_info['pop_names'])




# outfile=options.outfile
# if options.outfile_from_seqname:
#     name_without_parent_suffix="_".join(all_names[0][0].split("_")[:-1])
#     if outfile != "immediate_ancestry_results.txt":
#         outfile=name_without_parent_suffix+outfile+".txt"
#     else:
#         outfile=name_without_parent_suffix+".txt"



    
#for e,v in extra_info.items():
#    print(e,v)
    
#print(short_to_full)

#    print('seqs',sequences)
print('pops', extra_info['pop_names'])
#print('all_allele_frequencies',allele_frequencies)
#print('recombs',recombination_map)

#print(recombs)
#print(sequences)
#print(recombs)
# recombs_tmp=[]
# for recs in recombs:
#     recombs_tmp.append([0.1 for _ in recs])
# recombs=recombs_tmp
likelihoods=[generate_likelihood_from_data(allele_frequencies, 
                                           recombs, seq_system, 
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
elif options.type_of_analysis=='mcmc_search':
    popn=[full_to_short[n] for n in extra_info['pop_names']]
    ad=[]
    for likelihood in likelihoods:
        ad.append(mcmc_search(likelihood, popn, options.generations, short_to_full, N=options.mcmc_reps))
elif options.type_of_analysis=='evaluate_likelihoods':
    ad=[]
    for likelihood in likelihoods:
        ad.append(evaluate_list(likelihood, options.configs_to_test, options.generations, short_to_full))

res=''

print(ad)

for bd,outfile in zip(ad, options.outfiles):
    res=''
    for params,likel in bd:
        #print(params)
        #print(likel)
        full_names=list(map(str,params))
        short_names=[full_to_short[name] if len(name)>1 else name for name in full_names]
        res+=" ".join(short_names)+" "+str(likel)+"\n"
    print(res)
    with open(outfile   , "w") as f:
        f.write(res)


