import numpy


def transition_matrix(rho):
    '''
        input rho is the per-branch recombination probability, that is the probability that a recombination occurs on (only) one branch.
        input n is the number of generations back in time. 
    '''
    
    
    
    one_ancestry=numpy.array([
                        [1-2*rho,rho,rho/2.0,rho/2.0],
                        [rho,1-2*rho,rho/2.0,rho/2.0],
                        [rho/2.0,rho/2.0,1-2*rho,rho],
                        [rho/2.0,rho/2.0,rho,1-2*rho]])
    two_ancestries=numpy.kron(one_ancestry, one_ancestry)
    return two_ancestries

def generate_transition_matrix(rho_distances):
    def generator(index_of_sequence):
        return transition_matrix(rho_distances[index_of_sequence-1])
    return generator
    
def generate_emission_matrix(ancestral_allele_dictionary):
    '''
    The ancestral allele dictionary is a dictionary of the form {individual_number:[]}
    '''
    
    def subst(p,n):
        if n==1:
            return p
        return 1.0-p
    def generate(params):
        pops_mother=params[:4]
        pops_father=params[4:]
        def generator(index_of_sequence):
            ans=numpy.zeros((16,3))
            for i in range(4):   #the four ancestries of the mother
                for j in range(4): #the four ancestries of the father
                    for n1 in range(2): #the two alleles for one of the phases
                        for n2 in range(2): #the two alleles for the other phase
                            i1,i2= i*4+j, n1*2+n2
                            
                            #jumping over the unnecessary part.
                            if i2==2:
                                i2=1
                            if i2==3:
                                i2=2
                                
                            p1,p2= ancestral_allele_dictionary[pops_mother[i]][index_of_sequence], ancestral_allele_dictionary[pops_father[j]][index_of_sequence]
                            prob=0.5*subst(p1,n1)*subst(p2,n2)+0.5*subst(p1,n2)*subst(p2,n1)
                            ans[i1,i2]=prob
            return ans
        return generator
    
    return generate
                        
                        
                
            


if __name__=="__main__":
    alleles={"pop1":"ACGTT", "pop2":"CCCCC"}
    popsm=["pop1","pop2","pop1","pop1"]
    popsf=["pop2","pop2","pop2","pop2"]
    mu=0.01
    ad=generate_emission_matrix(alleles, mu)
    ad=ad(popsm+popsf)
    print(numpy.sum(ad(0),axis=0))
    print(numpy.sum(ad(1),axis=1))
    ad=generate_transition_matrix([0.1,0.2,0.3,0.4])
    print(numpy.sum(ad(0),axis=0))
    print(numpy.sum(ad(1),axis=1))
    
    