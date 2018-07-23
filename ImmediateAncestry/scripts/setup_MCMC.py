from numpy.random import choice
from summary import basic_summary
from meta_proposal import MetaProposal
from MCMC import basic_chain
from _operator import itemgetter
from ggrandparents_model import find_smallet_equivalence_class

def mcmc_search(likelihood, recombs, sequences, pops, generations, short_to_full, init=None, N=1000):
    print("#mcmc_serach")
#------------- Select of a random root to start
    print ("\nEscogemos nodo de inicio\n*********************************")
    if init is None:
        print("init = None")
        init=choice(pops,2**generations) 
        print ("init = ", init)
        #init=choice(pops,2)  # ['s' , 'e']
    init=''.join(map(str,init)) # se
    print ("*Nodo incial=", init)

#------------- 
	



    #-------------

    print("\n\tsummaries Llamda a basic_summary en summary.py\n\t*********************************")

    summaries=[basic_summary('posterior'), basic_summary('x')]
    summary_verbose_scheme={summary.name:(1,1) for summary in summaries}

    for key in summary_verbose_scheme:
        print ('\t* summary_verbose_scheme= ',key, ': ' , summary_verbose_scheme[key])


    print("\n\tproposal_function Llamda a MetaProposal en proposal.py\n\t*********************************")
    proposal_function=MetaProposal(pops=pops)
    print ("\t* proposal_function.Permutacion",proposal_function.props[0].pops)
    print ("\t* proposal_function.Sustitucion",proposal_function.props[1].pops)

#    print ("agfghj,kolkjmhngbfvd",sequences)

    print("\n\t _,_,s Llamda a basic_chain en MCMC.py\n\t*********************************")
    _,_, s=basic_chain(init, summaries, likelihood, proposal_function, post=None, N=N, 
                sample_verbose_scheme=summary_verbose_scheme, overall_thinning=1, i_start_from=0, 
                temperature=1.0)

#--------------
    res=list(s)
    print ("res=list(s)", res)
    #we know that the first column is iteration number, the second posterior, and the last is the configuration x.
    res_dic={find_smallet_equivalence_class(x):p for x,p in zip(res[2],res[1])}

    for key in res_dic:
        print ('* res_dic= ',key, ': ' , res_dic[key])



    xs, posteriors= map(list,list(zip(*list(res_dic.items()))))
    list_sorted=sorted(enumerate(posteriors), key=lambda elem: elem[1])
    i=list_sorted[-1][0]
    print("*xs[i] = ",xs[i])
    print ("*posterors[i] = ", posteriors[i])
    print("list_sorted", list_sorted)
    print("i=list_sorted[-1][0] = ", i)

    x_to_return=[]
    p_to_return=[]
    for i_tup in list_sorted[-10:][::-1]:
        i=i_tup[0]
        print("i=i_tup[0]",i)
        x_to_return.append(xs[i])
        print("x_to_return.append(xs[i]",x_to_return)
#        print ("XSSSSSSSSSSSSSSSSSSSSsss,",xs[i])
        p_to_return.append(posteriors[i])
        print("p_to_return.append(posteriors[i]", p_to_return)
    print("return list(zip(x_to_return, p_to_return)", x_to_return, p_to_return,"\n")
    return list(zip(x_to_return, p_to_return))
    
    
    
