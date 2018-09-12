def rescale(recombs, factor):
    new_recombs=[]
    for recomb_seq in recombs:
        new_recomb_seq=[]       
        for r in recomb_seq:
            new_recomb_seq.append(r*factor)
        new_recombs.append(new_recomb_seq)
    return new_recombs

