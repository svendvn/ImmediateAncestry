from numpy.random import choice, randint

def full_sequence(p):
    for r1,r2 in zip(p[::2],p[1::2]):
        if r1==r2:
            return False
    return True

def pure_sequence(p):
    return len(set(list(p)))==1

def simulate_two(seq, pops):
    n=len(seq)//2
    pop1,pop2=choice(pops,2, replace=False)
    new_seq=[pop1]*n+[pop2]*n
    return new_seq
    
    
    
def get_slice(parts, n):
    start=0
    end=n
    for i in parts:
        l=(end-start)//2
        start=start+i*l
        end=end-(1-i)*l
    return start, end

def sim_config(pops, complexity, generations=3):
    base_pop=choice(pops,1)[0]
    assert complexity<2**generations, 'The complexity is too high'
    sequence=[base_pop]*(2**generations)
    for i in range(complexity):
        candidate_sequence=sequence[:]
        n=len(candidate_sequence)
        parts=[]
        while True:
            if pure_sequence(candidate_sequence):
                break
            part=randint(0,2)
            n=n//2
            cand_candidate_sequence=candidate_sequence[part*n:(part+1)*n]
            if full_sequence(cand_candidate_sequence):
                part=1-part
                candidate_sequence=candidate_sequence[part*n:(part+1)*n]
            else:
                candidate_sequence=cand_candidate_sequence
            parts.append(part)
            
        n=len(sequence)
        #print parts
        start_slice, end_slice=get_slice(parts, len(sequence))
        #print start_slice, end_slice
        #print candidate_sequence
        #print 'sequence', sequence
        sequence[start_slice : end_slice]=simulate_two(candidate_sequence, pops)
    return sequence

if __name__=='__main__':
    pops=list('estv')
    print(sim_config(pops,7,4))
    