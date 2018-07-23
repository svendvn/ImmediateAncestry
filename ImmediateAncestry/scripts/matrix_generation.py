import numpy

from New_transition_emission_Matrix import transformation_EmissionMatrix
from New_transition_emission_Matrix import transformation_TransitionMatrix


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

def generate_transition_matrix(rho_distances, s, rho_infinity=False, max_numberof_recombs=-1, num_recombs=2):
#max_numberof_recombs == -1: OldTransMat; !=-1: NewTransMat
#num_recombs: Max num de recombs allowed. Default is 2 per chromosome.

    if rho_infinity:

       def inf_rho_generator(index_of_sequence):

            a=infinity_transition_matrix(s)

            if max_numberof_recombs == -1:
                return a
            else:
                return transformation_TransitionMatrix(a,num_recombs)

       return inf_rho_generator

    else:

       def generator(index_of_sequence):
          
          a=transition_matrix(rho_distances[index_of_sequence-1], s)

          if max_numberof_recombs == -1:
               return a

          else:

               return transformation_TransitionMatrix(a,num_recombs)

       return generator
    





def generate_emission_matrix(ancestral_allele_dictionary,s, max_numberof_recombs=-1, num_recombs=2 ):
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
    states=[]
    def generate(params):
        pops_mother=params[:m]
        pops_father=params[m:]
        
        #print("pops_mother,pops_father",pops_mother,pops_father)
        def generator(index_of_sequence):

            ans=numpy.zeros((M,4)) 
           # ans=numpy.zeros((M,5))
	 
            for i in range(m):   #the m ancestries of the mother
                for j in range(m): #the m ancestries of the father
                    p1,p2= ancestral_allele_dictionary[pops_mother[i]][index_of_sequence], ancestral_allele_dictionary[pops_father[j]][index_of_sequence]
                    i1=i*m+j
                    #print(pops_mother[i],pops_father[j])
                    estado=pops_mother[i]+" "+pops_father[j]
                    states.append(estado)
                    for n1 in range(2): #the two alleles and the missing for one of the phases
                        
                        for n2 in range(2): #the two alleles for the other phase
                            i2=n1*2+n2
                            #print (pops_mother[i],pops_father[j])
                            #jumping over the unnecessary part.
                            if i2==2:
                                continue
                            elif i2==3:
                                i2=2

                            prob= 0.5 * subst(p1,n1) * subst(p2,n2) + 0.5 * subst(p1,n2) * subst(p2,n1)
                            ans[i1,i2]=prob
                    #       states.append([pops_mother[i],pops_father[j]])
                    ans[i1,1]*=2
                    ans[i1,3]=1
                 #   print states[j]
                 #   ans[i1,4]=states[j]
         #   print(states[1])
         #   print (ans)

            if max_numberof_recombs == -1:
                return ans/2 #,states
            else:
                return transformation_EmissionMatrix(ans/2, int(num_recombs))
        return generator
    return generate
                               
            


if __name__=="__main__":
    alleles={"pop1":[0.2,0.1,0.1], "pop2":[0.9,0.9,0.9]}
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
    
    
