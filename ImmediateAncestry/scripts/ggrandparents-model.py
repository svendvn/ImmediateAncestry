from argparse import ArgumentParser
from numpy.random import choice
from ggrandparents_model import generate_likelihood_from_data
from simulate_data import simulate, simulate_recombs, simulate_allele_frequencies
from recombination import read_recombination_file
from construct_seqs import get_seqs
from brute_force_maximization import maximize_likelihood_exhaustive
from likelihood_evaluations import evaluate_list
from shortcut_names import read_shortcuts

usage="""This program generates a likelihood given the inputs: sequence, allele frequencies, and genetic, recombinational distances.

It then searches the spaces of possible ancestries of the great grandparents of the input sequence.

The data pipeline is used for testing purposes because the program contains ways to simulate data according to the model. There are basically four types of data to use:
1: no data at all, allele frequencies are simulated. It requires that you set the number of populations.
2: allele files of the ancestors are specified. 
3. The ancestors are simulated from the true_pops. It is possible to have more than one 

"""

parser = ArgumentParser(usage=usage)

#input/output arguments
parser.add_argument("--seq_files", type=str, nargs="+", default=['/home/svendvn/Dropbox/Bioinformatik/Zoomodel/data/prepared_data/seq12.txt',
                                                                 '/home/svendvn/Dropbox/Bioinformatik/Zoomodel/data/prepared_data/seq13.txt'], help="This is the string of a file of 0, 1 and 2s where 0 corresponds to two major allele, 1 corresponds to one of each and 2 is two minor alleles.")
parser.add_argument("--recomb_map", type=str, default=['/home/svendvn/Dropbox/Bioinformatik/Zoomodel/data/prepared_data/rho_12.txt',
                                                       '/home/svendvn/Dropbox/Bioinformatik/Zoomodel/data/prepared_data/rho_13.txt'], nargs="+", help="This is a list of files containing the genetic distances between SNPs. It is a tab separated line of numbers.")
parser.add_argument('--outfiles', type=str, nargs='+', default=["immediate_ancestry_results.txt"], help="This is the file in which the results are being stored.")
parser.add_argument("--allele_frequencies", type=str, default=['/home/svendvn/Dropbox/Bioinformatik/Zoomodel/data/prepared_data/chr12_freqs.txt',
                                                               '/home/svendvn/Dropbox/Bioinformatik/Zoomodel/data/prepared_data/chr13_freqs.txt'], nargs="+", help="This is a list of files containing the allele frequencies. They are in the format ")

#key arguments
parser.add_argument('--generations', type=int, default=3, help='the number of generations to go back. 3 is great grandparents.')

#technical arguments
parser.add_argument("--seq", type=str, default=[], nargs="+", help="If seq_files contains many lines of corresponding to different ancestors, this can be used to take out special rows. An empty list defaults to first and second row.")
parser.add_argument('--seq_indices', type=int, default=[0,1], nargs='+', help='the same as seq but here one can put the indices instead of the individualnames')
parser.add_argument('--ploidy_discrepancy', type=int, default=2, help='If the sequence files contain many different haploid sequences but you want diploid (unphased) data, set this to two and put two sequences in the --seq or --seq_indices option.')
parser.add_argument("--pops_to_search", type=str, nargs="+",default=[], help="This is a list of populations to search for the best solution. Defaults to all available population.")
parser.add_argument('--outfile_from_seqname', action='store_true', default=False, help="This will construct a filename from the names within the seq files and append the outfile-string and then '.txt'.")
parser.add_argument('--truncate_af', type=float, default=0.01, help='This will truncate the allele frequencies such that the lowest value is truncate_af and the highest value is 1-truncate_af')

#simulation arguments
parser.add_argument('--sequences_pipeline', type=int, nargs='+', default=[6,7], help='This is how to make the sequences to run the analysis on.')
parser.add_argument('--simulate_recombinations', default=False, action='store_true', help='this will simulate the recombination rate between SNPs. The upper rate is controlled by the argument recomb_rate')
parser.add_argument('--recomb_rate', type=float, default=0.1, help="The average recombination rate when simulated.")
parser.add_argument('--skewness', type=float, default=4, help="The skewness of the simulated recombination rates. If set to 0, there will be a constant recombination rate. Values below -1 are not meaningful")
parser.add_argument("--ancestors", type=int, default=[0,1,2,3,4,5,6,7], nargs="+", help="This list of numbers are the indices of the ancestors in the ancestor files which are chosen as ancestors.")
parser.add_argument("--ancestor_files", type=str, default=[], nargs="+", help="This is a list of files containing the ancestor alleles. Each file contains the alleles for a list of specimens for a certain segment of DNA \
                                                                                The file consists of the names of the specimen on the uneven lines and the alleles, space-separated, on the even lines.\
                                                                                If set to other than simulated, there will only be an effect if seqs is also simulated.")

parser.add_argument("--SNPs", type=int, default=0, help="this is the number of SNPs to simulate. If any input files are specified this will be ignored")
parser.add_argument("--reps", type=int, default=1, help="This is the number of independent segments of SNPs to draw if simulated")
parser.add_argument("--true_pops", type=str, nargs="+", default=[], help="If simulations take place, this will be the true immediate ancestors. It has to be of length 2**generations*ancestor_multiplier")
parser.add_argument('--no_pops', type=int, default=4, help='the number of populations.')
parser.add_argument('--ancestor_multiplier', type=int, default=1, help='the number of sets of ancestors are simulated. If more than one, the analysis will also be run more than once.')
parser.add_argument('--sequence_multiplier', type=int, default=1, help='the number of sequences that should be simulated from each set of ancestors.')
parser.add_argument('--population_names', type=str, nargs='+', default=[], help='When allele frequencies are simulated, the population names can not be read from anywhere, so to avoid default values (a,b,c,d, etc.) it can be specified here.')

parser.add_argument('--type_of_analysis', type=str, choices=['brute-force', 
                                                             'simulated_annealing', 
                                                             'evaluate_likelihoods'], 
                    default='evaluate_likelihoods',
                    help='chooses the type of analysis to run on the data. \
                          brute-force searches all possible combinations of gparents (based on the populations specified in allele_frequencies(that may be simulated)) \
                          simulated annealing is not implemented yet. \
                          evaluate_likelihoods evaluates the likelihoods specified in the list --configs_to_test')
parser.add_argument('--configs_to_test', type=str, nargs='+', default=['trivial_ellioti.txt'], help='if type of analysis is specified. If a string contains a dot, it will be read as filename')
parser.add_argument('--shortcut_names', type=str, default='shortcut_names.txt', help='file that short cuts long names for easier readability. It is of the form [full_name short_name\n,...]')

options = parser.parse_args()

if options.shortcut_names:
    full_to_short,short_to_full=read_shortcuts(options.shortcut_names)
else:
    identity={n:n for n in extra_info['pop_names']}
    full_to_short, short_to_full=identity

options.short_to_full=short_to_full

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

assert len(sequences)==len(options.outfiles)

assert set(full_to_short.keys())==set(extra_info['pop_names'])




# outfile=options.outfile
# if options.outfile_from_seqname:
#     name_without_parent_suffix="_".join(all_names[0][0].split("_")[:-1])
#     if outfile != "immediate_ancestry_results.txt":
#         outfile=name_without_parent_suffix+outfile+".txt"
#     else:
#         outfile=name_without_parent_suffix+".txt"



    

#    print('seqs',sequences)
print('pops', extra_info['pop_names'])
#print('all_allele_frequencies',allele_frequencies)
#print('recombs',recombination_map)

#print(recombs)
likelihoods=[generate_likelihood_from_data(allele_frequencies, recombs, seq_system, options.generations) for seq_system in sequences]
if options.type_of_analysis=='brute-force':
    ad=[]
    for likelihood in likelihoods:
        ad.append(maximize_likelihood_exhaustive(likelihood, pops, options.generations))
elif options.type_of_analysis=='simulated_annealing':
    assert False, 'not implemented'
elif options.type_of_analysis=='evaluate_likelihoods':
    ad=[]
    for likelihood in likelihoods:
        ad.append(evaluate_list(likelihood, options.configs_to_test, options.generations, short_to_full))

res=''



for bd,outfile in zip(ad, options.outfiles):
    res=''
    for params,likel in bd:
        full_names=list(map(str,params))
        short_names=[full_to_short[name] for name in full_names]
        res+=" ".join(short_names)+" "+str(likel)+"\n"
    print(res)
    with open(outfile   , "w") as f:
        f.write(res)


