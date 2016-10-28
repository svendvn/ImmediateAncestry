from argparse import ArgumentParser
from numpy.random import choice
from ggrandparents_model import generate_likelihood_from_data, maximize_likelihood_exhaustive
from simulate_data import simulate, simulate_recombs, simulate_allele_frequencies

usage="""This program generates a likelihood given the inputs: sequence, allele frequencies, and genetic, recombinational distances.

It then searches the spaces of possible ancestries of the great grandparents of the input sequence."""

parser = ArgumentParser(usage=usage)

parser.add_argument("--seq", type=str, nargs="+", default="simulated", help="This is the string of a file of 0, 1 and 2s where 0 corresponds to two major allele, 1 corresponds to one of each and 2 is two minor alleles.")
parser.add_argument("--alleles", type=str, default="simulated", nargs="+", help="This is a list of files containing the allele frequencies. The first line of each file contains the population name")
parser.add_argument("--recomb_map", type=str, default="simulated", help="This is a file containing the genetic distances between SNPs. It is a tab separated line of numbers.")
parser.add_argument("--no_pops", type=int, default=2, help="If not indicated elsewhere this will be the number of populations.")
parser.add_argument("--pops_to_search", type=str, nargs="+",default=[], help="This is a list of populations to search for the best solution. Defaults to all the populations from the supplied allele files.")
parser.add_argument("--SNPs", type=int, default=1000, help="this is the number of SNPs to simulate")
parser.add_argument("--reps", type=int, default=1, help="This is the number of independent segments of SNPs to draw if simulated")
parser.add_argument("--true_pops", type=str, nargs="+", default=[], help="If simulations take place, this will be the true immediate ancestors.")
parser.add_argument('--outfile', type=str, default="immediate_ancestry_results.txt", help="This is the file in which the results are being stored.")
parser.add_argument('--skewness', type=float, default=4, help="to be announced")
parser.add_argument('--recomb_rate', type=float, default=0.1, help="tba")

options = parser.parse_args()


setups=[]
pops=[]
if options.alleles!="simulated":
    no_pops=len(options.alleles)
    all_allele_frequencies=[]
    for n,r in enumerate(options.alleles):
        allele_frequencies={}
        lengths=[]
        with open(r, "r") as f:
            text=list(f.readlines())
            for i,j in zip(text[0::2], text[1::2]):
                name=i.rstrip()
                freqs=map(float, j.split())
                allele_frequencies[name]=freqs
                if n==0:
                    pops.append(name)
            lengths.append(len(freqs))
        all_allele_frequencies.append(allele_frequencies)
        setups.append(lengths)
if options.recomb_map != "simulated":
    lengths=[]
    recombs=[]
    for r in options.recomb_map:
        with open(options.recomb_map, "r") as f:
            recombs.append(map(float, f.readline().split()))
        lengths.append(len(recombs)+1)
if options.seq!= "simulated":
    seqs=[]
    lengths=[]
    for r in options.seq:
        with open(options.seq, "r") as f:
            seq=map(int, f.readline().split())
            lengths.append(len(seq))
        seqs.append(seq)
    setups.append(lengths)

res=""
print(setups)
if setups:
    setup=setups[0]
    reps=len(setup)
else:
    setup=[options.SNPs]*options.reps
    reps=options.reps

if options.alleles=="simulated":
    if options.pops_to_search:
        pops=set(pops_to_search)
    elif options.true_pops:
        pops=set(options.true_pops)
    else:
        pops=["pop"+str(i) for i in range(1,options.no_pops+1)]
    all_allele_frequencies= [simulate_allele_frequencies(pops, length) for length in setup]
if options.recomb_map=="simulated":
    recombs=[simulate_recombs(length,options.recomb_rate, options.skewness) for length in setup]
if options.seq=="simulated":
    if options.true_pops:
        tpops=options.true_pops
        assert len(tpops)==8, "true population not fully specified"
    else:
        tpops=[choice(pops) for _ in range(8)]
    seqs=[simulate(freq, tpops, recomb) for freq,recomb in zip(all_allele_frequencies, recombs)]
    

#print(recombs)


likelihood=generate_likelihood_from_data(all_allele_frequencies, recombs, seqs)
ad=maximize_likelihood_exhaustive(likelihood, pops)

if options.seq=="simulated":
    res=res+" ".join(map(str,tpops))+" "+ str(likelihood(tpops))+  "\n"
else:
    res=res+"#"+"\n"
print(res)

for params,likel in ad:
    res+=" ".join(map(str,params))+" "+str(likel)+"\n"

print(res)

with open(options.outfile, "w") as f:
    f.write(res)


