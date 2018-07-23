from numpy.random import random
from math import exp, log



            
def one_jump(x, post, temperature, posterior_function, proposal, pks={}):
    print ("#one_jump")
    
    newx,g1,g2,Jh,j1,j2=proposal(x,pks)
    print("newx,g1,g2,Jh,j1,j2 = ", newx,g1,g2,Jh,j1,j2)
    pks['g1']=g1
    pks['g2']=g2
    pks['Jh']=Jh
    pks['j1']=j1
    pks['j2']=j2
    
    post_new=posterior_function(newx,pks)
    print ("\t\tpost_new = ",post_new)

    pks['proposed_posterior']=post_new
    pks['proposed_x']=newx
    pks['old_x']=x
    pks['old_posterior']=post
  
    for key in pks:
        print ('\t\t pks InIcIaLiZaDo = ',key, ':' , pks[key])

#    print ("temperature"), temperature)



    #print 'Posterior Matrix:'
    #print (likelihood_old, prior_old)
    #print likelihood_new, prior_new
    
    #for key,val in pks.items():
    #    print key, '=', val
    
    if g2<=0 or j2<=0:
        logmhr=-float('inf')
        print ("g2<=0 & j2<=0 -> logmhr = ",logmhr) 
    else:
        logmhr=post_new-post+log(g2)+log(j2)-log(j1)-log(g1)+log(Jh)
        print ("g2>0 & j2>0 -> logmhr = ",logmhr) 

    if logmhr>100:
        mhr=float('inf')
        print("logmhr>100 -> mrh = ", mhr)
    else:
        mhr=exp(logmhr)
        print("logmhhr<= 100 -> mhr = ", mhr)
        
#    print ('MCMC.py, one_jump, post_new',post_new)
#    print ('MCMC.py, one_jump, post', post)
#    print ('MCMC.py, one_jump, post_new-post', post_new-post)
#    print ('MCMC.py, one_jump, exp(post_new-post)',exp(post_new-post)
#    print ('MCMC.py, one_jump, mhr', mhr)

    pks['mhr']=mhr
    u=random()
    pks['U']=u
    #proposal.adapt(mhr, u, post_new, post, temperature)
    
    if u<mhr:
        pks['posterior']=post_new
        pks['x']=newx
        print ("u<mhr -> post_new = ", post_new)
        print ("u<mhr -> newx = ", newx)
        return newx,post_new

    pks['posterior']=post
    pks['x']=x

    for key in pks:
        print ('\t\tpks AcTuAlIzAdO',key, ':' , pks[key])
    print ( "return x, post = ", x, post ,"\n")
    return x,post


def basic_chain(start_x, summaries, posterior_function, proposal, post=None, N=10000, 
                sample_verbose_scheme=None, overall_thinning=1, i_start_from=0, 
                temperature=1.0, 
                appending_result_file=None, appending_result_frequency=10):

    print ("#basic_chain(start_x, summaries, posterior_function, proposal, post=None, N=10000,                sample_verbose_scheme=None, overall_thinning=1, i_start_from=0,temperature=1.0,appending_result_file=None, appending_result_frequency=10)")
    x=start_x
    print ('\t* var x = ',x)
    if post is None:
        post=posterior_function(x)
        print ('\t* post = ',post)

    iteration_summary=[]
    #print 'random', random()
    count=0
    from_count=0

    print ("\tappending_result_file = ", appending_result_file)
    if appending_result_file is not None:
        with open(appending_result_file, 'w') as f:
            f.write(",".join(['iteration'] + [s.name for s in summaries])+'\n')

        
    for i in range(i_start_from,i_start_from+N):
        print ("\t*i_start_from = ", i_start_from)
        print ("\ti_start_from+N = ", i_start_from+N)

        print ("\n\t\tLlamada one_jump en MCMC.py \n\t\t*********************************")
        proposal_knowledge_scraper={}
        new_x,new_post=one_jump(x, post, temperature, posterior_function, proposal, proposal_knowledge_scraper)
        print ("\tnew_x", new_x)
        print ("\tnew_post", new_post)



        if overall_thinning!=0 and i%overall_thinning==0:
            print("overall_thinning!=0 and i%overall_thinning==0")
            iteration_summary.append(_calc_and_print_summaries(sample_verbose_scheme,
                                                               summaries,
                                                               iteration_number=i,**proposal_knowledge_scraper))
            print ("\toveral_thinning = ", overall_thinning)
            print ("iteration_summary = ", iteration_summary)



            print ("appending_result_file = ", appending_result_file)
            if appending_result_file is not None:
                count+=1
                print ("appending_result_file = None -> count =  ", count)
                if count % appending_result_frequency==0:
                    with open(appending_result_file, 'a') as f:
                        for n,params in enumerate(iteration_summary[from_count:]):
                            f.write(",".join(map(str, params))+'\n')
                    from_count=count
                    print ("count % appending_result_frequency==0 -> form_count =  ", form_count)

        x=new_x
        post=new_post
    
    print("return x=new_x , post=new_post= , iteration_summay",x, post, *iteration_summary,"\n")

#    print("* iteration_summary = ",iteration_summary)
    return x, post, zip(*iteration_summary)
        
        
        
def _calc_and_print_summaries(sample_verbose_scheme,summaries,**kwargs):
    print ("#_calc_and_print_summaries(sample_verbose_scheme,summaries,**kwargs)")

    iteration=kwargs['iteration_number']
    print ("iteration", iteration)
    res=[iteration]
    print ("res = ", res)
    for s in summaries:
        print("s = " , s)
        save_num,print_num=sample_verbose_scheme.get(s.name, (0,0))
        print ("save_num,print_num = ", save_num,print_num)
        save_bool = (save_num!=0) and (iteration % save_num==0) 
        print ("save_bool", save_bool)
        print_bool = (print_num!=0) and (iteration % print_num==0)
        print("print_bool",print_bool)
        if save_bool or print_bool:
            print ("if save_bool or print_bool")
            val=s(**kwargs)
            print ("val = ", val)
            if print_bool:
                print ("if print_bool")               
                print(str(iteration)+'. '+ s.pretty_print(val))
            if save_bool:
                print ("if save_bool")
                res.append(val)
                print ("res.append(val) = ", res)
            else:
                res.append(None)
        else:
            res.append(None)
    print ("return res ", res, "\n")
    return res
    
