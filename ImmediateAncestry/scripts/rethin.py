def rethin_sequences(sequences, thin_coef=10):
    new_sequences=[]
    for sequence_system in sequences:
        chroms=[]
        for sequence in sequence_system:
            chroms.extend(rethin_seq_chrom(sequence, thin_coef))
        new_sequences.append(chroms)
    return new_sequences
            
def rethin_seq_chrom(sequence, thin_coef=10):
    new_chroms=[]
    for i in range(thin_coef):
        new_chroms.append(sequence[i::thin_coef])
    return new_chroms
    
def rethin_recombination(recombs, thin_coef=10):
    new_recombs=[]
    for rec_seq in recombs:
        new_recombs.extend(rethin_recomb_chrom(rec_seq, thin_coef))
    return new_recombs

def rethin_recomb_chrom(rec_seq, thin_coef):
    new_chroms=[]
    for i in range(thin_coef):
        j=i
        chrom=[]
        while j<len(rec_seq)-thin_coef+1:
            chrom.append(sum(rec_seq[j:j+thin_coef]))
            j+=thin_coef
        new_chroms.append(chrom)
    return new_chroms

def rethin_allele_frequencies(allele_frequencies, thin_coef=10):
    allele_freq_chroms=[]
    for chrom_dic in allele_frequencies:
        allele_freq_chroms.extend(rethin_allele_chrom(chrom_dic, thin_coef))
    return allele_freq_chroms

def rethin_allele_chrom(chrom_dic, thin_coef=10):
    chroms=[]
    for i in range(thin_coef):
        new_dic={}
        for pop,freqs in chrom_dic.items():
            new_dic[pop]=freqs[i::thin_coef]
        chroms.append(new_dic)
    return chroms

if __name__=='__main__':
    seq=[[[1,2,3,4,5,6,7,8,9],[0,1,2,3,4,5,6,7]]]
    thin_coef=3
    rec=[[0.1,0.3,0.2,0.1,0.05,0.15,0.22,0.2],[0.0,0.1,0.02,0.04,0.1,0.3,0.01]]
    allele_freqs=[
        {'t':[0.11,0.31,0.21,0.11,0.051,0.151,0.221,0.21,0.1]
            },
        {'t':[0.01,0.11,0.021,0.041,0.11,0.31,0.011,0.8]}]

    print(rethin_sequences(seq,thin_coef))
    print(rethin_allele_frequencies(allele_freqs, thin_coef))
    print(rethin_recombination(rec, thin_coef))
     