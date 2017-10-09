from ggrandparents_model import find_smallet_equivalence_class
from itertools import product

def maximize_likelihood_exhaustive(likelihood, pops_to_choose_from, generations=3):
    combinations=[]
    likelihoods=[]
    counter=0
    n=2**generations
    print('pops_to_choose_from', pops_to_choose_from)
    for itera in product(*([pops_to_choose_from]*n)):
        if find_smallet_equivalence_class(itera) in combinations:
            print(str(itera), "skipped")
            continue
        counter+=1
        combinations.append(itera)
        #likelihoods.append(counter)
        print(itera)
        likelihoods.append(likelihood(list(itera)))
    sorted_indexes=[i[0] for i in sorted(enumerate(likelihoods), key=lambda x: x[1])]
    res_dict=[]
    for i in range(len(likelihoods)):#print everything.
        index=sorted_indexes[-i-1]
        res_dict.append((combinations[index],likelihoods[index]))
    print("Looped over ",counter)
    return res_dict