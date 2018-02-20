import numpy

REORDER={0:0,1:1,2:2,4:3,5:4,8:5}

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

def generate_transition_matrix(rho_distances, s):
    def generator(index_of_sequence):
        return transition_matrix(rho_distances[index_of_sequence-1], s)
    return generator
    
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
    def generate(params):
        pops_mother=params[:m]
        pops_father=params[m:]
        def generator(index_of_sequence):
            ans=numpy.zeros((M,6))
            for i in range(m):   #the m ancestries of the mother
                for j in range(m): #the m ancestries of the father
                    for n1 in range(3): #the two alleles and the missing for one of the phases
                        for n2 in range(3): #the two alleles for the other phase
                            i1,i2= i*m+j, n1*3+n2
                            
                            #jumping over the unnecessary part.
                            if i2 not in REORDER:
                                continue
                            else:
                                i2=REORDER[i2]
                                
                            p1,p2= ancestral_allele_dictionary[pops_mother[i]][index_of_sequence], ancestral_allele_dictionary[pops_father[j]][index_of_sequence]
                            prob= 0.5*subst(p1,n1)*subst(p2,n2)+0.5*subst(p1,n2)*subst(p2,n1)
                            ans[i1,i2]=prob
                            
                                
            return ans
        return generator
    
    return generate
                        
                        
                
            


if __name__=="__main__":
    alleles={"pop1":[0.1,0.1,0.1], "pop2":[0.9,0.9,0.9]}
    popsm=["pop1","pop2","pop1","pop1"]
    popsf=["pop2","pop2","pop1","pop2"]
    mu=0.01
    ad=generate_emission_matrix(alleles,3)
    ad=ad(popsm+popsf)
    print(ad(0))
    print(numpy.sum(ad(0),axis=0))
    print(numpy.sum(ad(1),axis=1))
    ad=generate_transition_matrix([0.1,0.0,0.0,0.0],3)
    print(ad(1))
    print(numpy.sum(ad(0),axis=0))
    print(numpy.sum(ad(1),axis=1))
    
    #we expect to get fff,ffm,fmf,fmm,mff,mfm,mmf,mmm
    
    