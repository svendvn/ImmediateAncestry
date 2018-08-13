def get_val(listi):
    return "".join(listi)

def get_best_of_size(listi, size):
    assert len(listi)%size==0, "wrong list size"
    assert size%2==0, "wrong size"
    band=listi
    h= size//2
    bval= get_val(listi)
    for i in range(len(listi)//size):
        cand=band[:(size*i)]+\
                band[(size*i+h):(size*(i+1))]+\
                band[(size*i):(size*i+h)]+\
                band[(size*(i+1)):]
        
        cval=get_val(cand)
        if cval<bval:
            band=cand
            bval=cval
    return bval



def find_smallest_equivalence_class(params):
    size=2
    best_combi=params
    while size<=len(params):
        best_combi=get_best_of_size(best_combi, size)
        size*=2
    return best_combi

if __name__ == '__main__':
    #call_helper()
    #print(call_forward())
    #lik=test_model_likelihood(8)
    #print(maximize_likelihood_exhaustive(lik,["pop1","pop2"]))
    print(find_smallest_equivalence_class('eevt'))
    print(find_smallest_equivalence_class(list('eevt')))