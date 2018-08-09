
def bin_sequences(sequences, snps_per_bin):
    bin_maps=[]
    seq_system=sequences[0]
    for sequence in seq_system:
        bin_map=bin_uniformly(sequence, snps_per_bin)
        bin_maps.append(bin_map)
    return bin_maps

def bin_uniformly(sequence, snps_per_bin=100):
    n=len(sequence)
    bin_map=[]
    no_blocks=n//snps_per_bin
    if no_blocks>1:
        for i in range(no_blocks):
            bin_map.append(list(range((i*n)//no_blocks, ((i+1)*n)//no_blocks)))
    else:
        bin_map=list(range(n))
    return bin_map

def collapse_recombination_map(recombs, bin_maps):
    new_recombs=[]
    for recomb,bin_map in zip(recombs, bin_maps):
        new_recomb=[]
        for n,l in enumerate(bin_map):
            m=len(l)
            m2=m//2
            if n>0:
                new_recomb[-1]+=sum((recomb[i] for i in l[:m2]))
            if n<len(bin_map)-1:
                new_recomb.append(sum((recomb[i] for i in l[m2:])))
        new_recombs.append(new_recomb)
    return new_recombs
            

if __name__=='__main__':
    sequence=[0]*109
    recombs=[1]*109
    print(bin_uniformly(sequence,10))
    nr=collapse_recombination_map([recombs], [bin_uniformly(sequence, 10)])
    print(nr)