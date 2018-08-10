import numpy

REORDER={0:0,1:1,3:2}

def transition_matrix(rho, s):
    '''
        input rho is the per-branch recombination probability, that is the probability that a recombination occurs on (only) one branch.
        input s is the number of generations back in time where great grandparents are 3. 
    '''
    
    dig_element=1-(s-1)*rho
    elements=[rho/(2**k) for k in range(s-1)]
    one_ancestry=numpy.identity(2**(s-1))*dig_element
    for i in range(2**(s-1)):
        for j in range(i+1, 2**(s-1)):
            c=1
            while i//(2**c)!=j//(2**c):
                c+=1
            one_ancestry[i,j]=elements[c-1]
            one_ancestry[j,i]=elements[c-1]
    two_ancestries=numpy.kron(one_ancestry, one_ancestry)
    return two_ancestries

def infinity_transition_matrix(s):
    M=(2**(s-1))**2
    return numpy.ones((M,M))*(1.0/M)

def initial_matrix(generations):
    M=(2**(generations-1))**2
    M=[1.0/M]*M
    return M

def generate_transition_matrix(rho_distances, s, rho_infinity=False):
    if rho_infinity:
        def inf_rho_generator(index_of_sequence):
            return infinity_transition_matrix(s)
        return inf_rho_generator
    def generator(index_of_sequence):
        return transition_matrix(rho_distances[index_of_sequence-1], s)
    return generator

def calculate_before_dic(bin_map, ancestral_allele_frequencies, sequence, pops):
    avg_af=numpy.mean(numpy.array(list(ancestral_allele_frequencies.values())),axis=0)
    probabilities=[] #of the form [{'ee':0.331, 'et':0.1234,..., 'vv':0.0013}, {...},...,{...}].
    pop_combinations=[tuple(sorted([p1,p2])) for n,p1 in enumerate(pops) for m,p2 in enumerate(pops) if n>=m]
    for snp_list in bin_map:
        bin_probs={pop_comb:0 for pop_comb in pop_combinations}
        for i in snp_list:
            ps= {pop:ancestral_allele_frequencies[pop][i] for pop in pops}
            for pop_comb in pop_combinations:
                bin_probs[pop_comb]+=numpy.log(calc_prob(ps[pop_comb[0]], 
                                                         ps[pop_comb[1]], 
                                                         sequence[i]))
                bin_probs[pop_comb]-=numpy.log(calc_prob(avg_af[i],
                                                         avg_af[i],
                                                         sequence[i]))
        probabilities.append(bin_probs)
    #print(len(sequence))
    #print(probabilities)
    return probabilities


def subst(p,n):
    if n==1:
        return p
    return 1.0-p

def calc_prob(prob1,prob2,observed):
    if observed==3:
        return 1
    elif observed==2:
        n1=n2=1
    elif observed==1:
        n1=1
        n2=0
    elif observed==0:
        n1=n2=0
    prob= 0.5 * subst(prob1,n1) * subst(prob2,n2) + 0.5 * subst(prob1,n2) * subst(prob2,n1)
    if observed==1:
        return 2*prob
    return prob
                
        
    
def generate_emission_matrix_binned(before_dic,s):
    '''
    The ancestral allele dictionary is a dictionary of the form {individual_number:[]}
    s is the number of generations back such that s=3 is great grandparents
    '''
    
    

    m=2**(s-1)
    M=m**2
    def generate(params, log=False):
        pops_mother=params[:m]
        pops_father=params[m:]
        def generator(index_of_sequence):
            ans=numpy.zeros((M,1))
            for i,pop_mom in enumerate(pops_mother):   #the m ancestries of the mother
                for j,pop_dad in enumerate(pops_father): #the m ancestries of the father
                    i1=i*m+j
                    if log:
                        ans[i1,0]=before_dic[index_of_sequence][tuple(sorted([pop_dad,pop_mom]))]
                    else:
                        ans[i1,0]=numpy.exp(before_dic[index_of_sequence][tuple(sorted([pop_dad,pop_mom]))])
            return ans
        return generator
    
    return generate


def generate_emission_matrix(ancestral_allele_dictionary,s):
    '''
    The ancestral allele dictionary is a dictionary of the form {individual_number:[]}
    s is the number of generations back such that s=3 is great grandparents
    '''
    
    
    def subst(p,n):
        if n==1:
            return p
        elif n==2:
            return 0.5
        return 1.0-p
    m=2**(s-1)
    M=m**2
    def generate(params, log=False):
        pops_mother=params[:m]
        pops_father=params[m:]
        def generator(index_of_sequence):
            ans=numpy.zeros((M,4))
            for i in range(m):   #the m ancestries of the mother
                for j in range(m): #the m ancestries of the father
                    p1,p2= ancestral_allele_dictionary[pops_mother[i]][index_of_sequence], ancestral_allele_dictionary[pops_father[j]][index_of_sequence]
                    i1=i*m+j
                    for n1 in range(2): #the two alleles and the missing for one of the phases
                        for n2 in range(2): #the two alleles for the other phase
                            i2=n1*2+n2
                            
                            #jumping over the unnecessary part.
                            if i2==2:
                                continue
                            elif i2==3:
                                i2=2

                            prob= 0.5 * subst(p1,n1) * subst(p2,n2) + 0.5 * subst(p1,n2) * subst(p2,n1)
                            ans[i1,i2]=prob
                    ans[i1,1]*=2
                    ans[i1,3]=1
                            
            if log:
                ans=numpy.log(ans)                    
            return ans
        return generator
    
    return generate

                        
                        
                
            


if __name__=="__main__":
    alleles={"pop1":[0.0384], "pop2":[0.4722]}
    popsm=["pop2","pop2","pop2","pop2"]
    popsf=["pop2","pop2","pop2","pop2"]
    mu=0.01
    ad=generate_emission_matrix(alleles,3)
    ad=ad(popsm+popsf)
    print(ad(0))
    print(numpy.sum(ad(0),axis=0))
    print(numpy.sum(ad(0),axis=1))
    ad=generate_transition_matrix([0.1,0.0,0.0,0.0],3)
    print(ad(1))
    #print(numpy.sum(ad(0),axis=0))
    #print(numpy.sum(ad(1),axis=1))
    
    #we expect to get fff,ffm,fmf,fmm,mff,mfm,mmf,mmm
    
    sequence=[0]*9
    allele_freqs={'a':[0.1]*3+[0.2]*6,'b':5*[0.01]+2*[0.99]+2*[0.01]}
    bin_map=[list(range(4)), list(range(4,9))]
    generations=2
    params=['a','b','a','b']
    a=calculate_before_dic(bin_map, allele_freqs, sequence, pops=set(allele_freqs.keys()))
    print(bin_map)
    print(a)
    ad=generate_emission_matrix_binned(a,generations)
    ad_params=ad(params)
    print(ad_params(0))
    print(ad_params(1))
    