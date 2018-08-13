from numpy.random import choice
from summary import basic_summary
from meta_proposal import MetaProposal
from MCMC import basic_chain
from _operator import itemgetter
from configuration_classes import find_smallest_equivalence_class

def mcmc_search(likelihood, pops, generations, short_to_full, init=None, N=1000):
    '''
    This function wraps the MCMC call, and takes the MCMC output and finds the maximum posterior value.
    '''
    if init is None:
        init=choice(pops,2**generations)
        
    init=''.join(map(str,init))
        
    summaries=[basic_summary('posterior'), basic_summary('x')]
    summary_verbose_scheme={summary.name:(1,1) for summary in summaries}
    
    proposal_function=MetaProposal(pops=pops)
    
    _,_, s=basic_chain(init, summaries, likelihood, proposal_function, post=None, N=N, 
                sample_verbose_scheme=summary_verbose_scheme, overall_thinning=1, i_start_from=0, 
                temperature=1.0)
    res=list(s)
    #print(res)
    #we know that the first column is iteration number, the second posterior, and the last is the configuration x.
    res_dic={find_smallest_equivalence_class(x):p for x,p in zip(res[2],res[1])}
    xs, posteriors= map(list,list(zip(*list(res_dic.items()))))
    list_sorted=sorted(enumerate(posteriors), key=lambda elem: elem[1])
    i=list_sorted[-1][0]
    #print(xs[i], posteriors[i])
    x_to_return=[]
    p_to_return=[]
    for i_tup in list_sorted[-10:][::-1]:
        i=i_tup[0]
        x_to_return.append(xs[i])
        p_to_return.append(posteriors[i])

    return list(zip(x_to_return, p_to_return))
    
    
    