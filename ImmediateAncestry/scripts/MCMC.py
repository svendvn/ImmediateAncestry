from numpy.random import random
from math import exp, log



            
def one_jump(x, post, temperature, posterior_function, proposal, pks={}):
    
    newx,g1,g2,Jh,j1,j2=proposal(x,pks)
    pks['g1']=g1
    pks['g2']=g2
    pks['Jh']=Jh
    pks['j1']=j1
    pks['j2']=j2
    
    post_new=posterior_function(newx,pks)
    pks['proposed_posterior']=post_new
    pks['proposed_x']=newx
    pks['old_x']=x
    pks['old_posterior']=post
    
    #print temperature

    #print 'Posterior Matrix:'
    #print likelihood_old, prior_old
    #print likelihood_new, prior_new
    
    #for key,val in pks.items():
    #    print key, '=', val
    
    if g2<=0 or j2<=0:
        logmhr=-float('inf')
    else:
        logmhr=post_new-post+log(g2)+log(j2)-log(j1)-log(g1)+log(Jh)
    if logmhr>100:
        mhr=float('inf')
    else:
        mhr=exp(logmhr)
        
    #print post_new, post, post_new-post, exp(post_new-post), mhr
    pks['mhr']=mhr
    
    u=random()
    pks['U']=u
    #proposal.adapt(mhr, u, post_new, post, temperature)
    
    if u<mhr:
        pks['posterior']=post_new
        pks['x']=newx
        return newx,post_new
    pks['posterior']=post
    pks['x']=x
    return x,post


def basic_chain(start_x, summaries, posterior_function, proposal, post=None, N=10000, 
                sample_verbose_scheme=None, overall_thinning=1, i_start_from=0, 
                temperature=1.0, 
                appending_result_file=None, appending_result_frequency=10):
    
    x=start_x
    if post is None:
        post=posterior_function(x)
    
    iteration_summary=[]
    #print 'random', random()
    count=0
    from_count=0
    
    if appending_result_file is not None:
        with open(appending_result_file, 'w') as f:
            f.write(",".join(['iteration'] + [s.name for s in summaries])+'\n')
        
    for i in range(i_start_from,i_start_from+N):
        proposal_knowledge_scraper={}
        new_x,new_post=one_jump(x, post, temperature, posterior_function, proposal, proposal_knowledge_scraper)
        if overall_thinning!=0 and i%overall_thinning==0:
            iteration_summary.append(_calc_and_print_summaries(sample_verbose_scheme,
                                                               summaries,
                                                               iteration_number=i,**proposal_knowledge_scraper))

            if appending_result_file is not None:
                count+=1
                if count % appending_result_frequency==0:
                    with open(appending_result_file, 'a') as f:
                        for n,params in enumerate(iteration_summary[from_count:]):
                            f.write(",".join(map(str, params))+'\n')
                    from_count=count
        x=new_x
        post=new_post
    
    #print iteration_summary
    #print zip(*iteration_summary)
    return x, post, zip(*iteration_summary)
        
        
        
def _calc_and_print_summaries(sample_verbose_scheme,summaries,**kwargs):
    iteration=kwargs['iteration_number']
    res=[iteration]
    for s in summaries:
        save_num,print_num=sample_verbose_scheme.get(s.name, (0,0))
        save_bool = (save_num!=0) and (iteration % save_num==0) 
        print_bool = (print_num!=0) and (iteration % print_num==0)
        if save_bool or print_bool:
            val=s(**kwargs)
            if print_bool:
                print(str(iteration)+'. '+ s.pretty_print(val))
            if save_bool:
                res.append(val)
            else:
                res.append(None)
        else:
            res.append(None)
    return res
    