import numpy

max_recomb_rate=0.1

def simulate_all(length, pops=(["pops1","pops1","pops2","pops1"]+["pops2"]*4)):
    print(pops)
    allele_frequencies=simulate_allele_frequencies(set(pops), length)
    recombs=simulate_recombs(length-1)
    return simulate(allele_frequencies, pops, recombs)
    
def simulate_allele_frequencies(set_of_populations, length):
    res={}
    for j in set_of_populations:
        add=[]
        for _ in range(length):
            add.append(numpy.random.random())
        res[j]=add
    return res

def simulate_recombs(length,avg_recomb_rate, skewness=4):
    max_recomb_rate=(skewness+1)*avg_recomb_rate
    return [(numpy.random.random())**skewness*max_recomb_rate for _ in range(length)]
        

def simulate(recombs, allele_frequencies=None, pops=None, ancestors="simulated"):
    """
    Simulates one sequence in the end containing 0,1 and 2s. pops is the true parameters, that is the ancestry of the four files.
    """
    
    if ancestors=="simulated":
        assert pops is not None, "populations to simulate from not specified."
        assert allele_frequencies is not None, "allele frequencies to simulate from not specified."
        seqs=sim_ancestries(allele_frequencies, pops)
    
    else:
        seqs=ancestors
    
    #print(recombs)
    #print(ancestors)
    #print(seqs)
    while len(seqs)>2:
        print(len(seqs))
        seqs=sim_next_gen(seqs, recombs)
        print(len(seqs))
        #print(seqs)
    
    return [i+j for i,j in zip(seqs[0],seqs[1])]
    #write_down(res_file, seqs)
    
    
def sim_next_gen(seqs, recombs):
    
    assert len(seqs)%2==0, "unequal number of parents"
    res=[]
    
    #print(seqs)
    
    for m,(father,mother) in enumerate(zip(seqs[0::2], seqs[1::2])):
        child=[]
        inherit_from=int(numpy.random.random()<0.5) # random to inherit from both of them.
        print(len(recombs), len(father),len(mother))
        if len(recombs)<(len(father)-1):
            print("m",m)
        for n,fm_alleles in enumerate(zip(father,mother)):
            child.append(fm_alleles[inherit_from])
            if n<(len(father)-1):
                if numpy.random.random()<recombs[n]: #make recombination
                    inherit_from={0:1,1:0}[inherit_from]
        res.append(child)
        
    return res
    
def sim_ancestries(allele_frequencies, pops):
    res=[]
    for pop in pops:
        add=[]
        for prob in allele_frequencies[pop]:
            add.append(int(numpy.random.random()<prob))
        res.append(add)
    return res
    
def write_down(file_name, seq):
    with open(file_name, "w") as f:
        for i in seq:
            f.write(str(i)+" ")

if __name__=="__main__":
    allele_fs={"pop1":[0.1,0.2,0.3,0.4,0.04,0.03,0.1], "pop2":[0.8,0.5,0.99,0.4,0.3,0.99,0.7]}
    res_file="geschaft.txt"
    pops=["pop1","pop1","pop1","pop1","pop2","pop2","pop2","pop2"]
    recombs=[0.01,0.01,0.01,0.99,0.01,0.01]
    #simulate(allele_fs, res_file, pops, recombs)
    write_down(res_file,simulate_all(10000, pops))
    