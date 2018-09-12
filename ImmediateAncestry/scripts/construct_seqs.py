from simulate_data import simulate_allele_frequencies, sim_ancestries, simulate_recombs, simulate, zero_one_ize
from numpy.random import choice
from rethin import rethin_sequences, rethin_allele_frequencies, rethin_recombination
from snp_binner import collapse_recombination_map, bin_sequences
from rescale_recombinations import rescale


def rethin_wrapper(statistic, options, extra_info):
    sequences, allele_frequencies, recombs=statistic
    new_sequences=rethin_sequences(sequences,options.thin_coef)
    new_allele_frequencies=rethin_allele_frequencies(allele_frequencies, options.thin_coef)
    new_recombs=rethin_recombination(recombs, options.thin_coef)
    return new_sequences, new_allele_frequencies, new_recombs

def bin_wrapper(statistic, options, extra_info):
    sequences, allele_frequencies, recombs=statistic
    bin_maps=bin_sequences(sequences, options.bin_window_size)
    new_recombs=collapse_recombination_map(recombs,bin_maps)
    return sequences, allele_frequencies, new_recombs, bin_maps

def rescale_wrapper(statistic, options, extra_info):
    if len(statistic)==4:
        sequences, allele_frequencies, recombs, bin_maps=statistic
    else:
        sequences, allele_frequencies, recombs = statistic
    new_recombs=rescale(recombs, options.rescale_recombinations_factor)
    if len(statistic)==4:
        return sequences, allele_frequencies, new_recombs, bin_maps
    else:
        return sequences, allele_frequencies, new_recombs
    
    


def simulate_allele_frequencies_wrapper(statistic, options, extra_info):
    if 'setups' not in extra_info:
        extra_info['setups']=[options.SNPs]*options.reps
    allele_frequencies=[] #about to be filled with dictionaries
    for length in extra_info['setups']: 
        this_segment_allele_frequencies=simulate_allele_frequencies(statistic, length)
        allele_frequencies.append(this_segment_allele_frequencies)
    extra_info['allele_frequencies']=allele_frequencies
    return allele_frequencies

def simulate_ancestors_wrapper(all_allele_frequencies, options, extra_info):
    res=[]
    if options.true_pops is None:
        r=[]
        for _ in  range(2**options.generations):
            print(extra_info['pop_names'])
            i=choice(len(extra_info['pop_names']),1)
            a=extra_info['pop_names'][i]
            r.append(a)
        options.true_pops=[r]
    for pops in options.true_pops:
        tmp=[]
        for allele_frequencies in all_allele_frequencies:
            tmp.append(sim_ancestries(allele_frequencies, [options.short_to_full[pop] for pop in pops]))
        res.append(tmp)
    return res

def load_recombinations_wrapper(statistic, options, extra_info):
    if 'recombs' in extra_info:
        return (statistic, extra_info['recombs'])
    elif options.simulate_recombination:
        res=[]
        for length in extra_info['setups']:
            res.append(simulate_recombs(length,options.sim_recomb_rate, skewness=options.skewness))
        extra_info['recombs']=res
        return (statistic,res)
    else:
        assert options.simulate_recombination, 'It seems that recombination map should neither be loaded nor simulated. Therefore, error.'
            
def simulate_generations_wrapper(all_ancestors_and_all_recombs, options, extra_info):
    res=[]
    true_vals=[]
    all_ancestors, all_recombs=all_ancestors_and_all_recombs
    for ancestor_tuple in all_ancestors:
        for _ in range(options.sequence_multiplier):
            tmp=[]
            for ancestor_piece, recomb_piece in zip(ancestor_tuple,all_recombs):
                tmp.append(simulate(recombs=recomb_piece, ancestors=ancestor_piece))
            res.append(tmp)
            true_vals.append(ancestor_tuple)
    extra_info['true_vals']=true_vals
    return res



def package_the_deal_wrapper(statistic, options, extra_info):
    if 'allele_frequencies' not in extra_info:
        _=read_allele_frequencies(options.allele_frequencies, options.truncate_af, extra_info)
    return (statistic, extra_info['allele_frequencies'], extra_info['recombs'])


   
    

            
    

def get_indices_of_pops(chosen_pops, all_names, short_to_full):
    res=[]
    for chosen_pop in chosen_pops:
        tmp=[]
        for index,name in enumerate(all_names):
            if name.startswith(short_to_full[chosen_pop]):
                tmp.append(index)
        res.append(tmp)
    return res

def randomize_inds(list_of_lists):
    res=[]
    for indices in list_of_lists:
        res.append(choice(indices))
    return res

def choose_ancestors_wrapper(ancestors, options, extra_info):
    '''
    This should either
    a) extract fully specified ancestors in the options.ancestors or options.ancestor_indices.This is chosen if the length of those lists 
    matches what it should and a warning is sent out if ##NOT IMPLEMENTED!
    b) simulate random ancestors from the shortcut names object called options.short_to_full
    and the true populations called true_pops.
    '''
    if options.ancestors or options.ancestor_indices:
        assert False, 'Not implemented case sought after'
    all_chosen_ancestors=[]
    #print('short_to_full', options.short_to_full)
    
    for true_pop_config in options.true_pops:
        chosen_ancestors=[]
        #print('true_pop_config', true_pop_config)
        inds= get_indices_of_pops(true_pop_config, extra_info['names'], options.short_to_full)
        randomized_inds=randomize_inds(inds)
        #print('rinds', randomized_inds)
        for ancestor_pieces in ancestors:
            chosen_ancestors_one_piece=[]
            for ind in randomized_inds:
                chosen_ancestors_one_piece.append(zero_one_ize(ancestor_pieces[ind]))
            chosen_ancestors.append(chosen_ancestors_one_piece)
        all_chosen_ancestors.append(chosen_ancestors)
    #print(all_chosen_ancestors)
    return all_chosen_ancestors
            
            
        
    
    assert False, 'stopped implementing'

transitions={
             (1,3): simulate_allele_frequencies_wrapper,
             #(2,3): calculate_allele_frequencies_wrapper,
             (3,4): simulate_ancestors_wrapper,
             (2,4): choose_ancestors_wrapper,
             (4,5): load_recombinations_wrapper,
             (5,6): simulate_generations_wrapper,
             #(3,6): simualte_mixed_individual_wrapper,
             #(2,6): choose_seqs_wrapper
             (6,7): package_the_deal_wrapper,
             (7,8): rethin_wrapper,
             (7,9): bin_wrapper,
             (8,9): bin_wrapper,
             (7,10): rescale_wrapper,
             (8,10): rescale_wrapper,
             (9,10): rescale_wrapper
             }



def get_seqs(options, extra_info={}):
    pipeline=options.sequences_pipeline
    
    statistic=get_input(pipeline[0], options, extra_info)
    previous_stage=pipeline[0]
    
    for stage in pipeline[1:]:
        trans=(previous_stage, stage)
        assert trans in transitions, 'The transition is illegal or not implemented yet'
        statistic= transitions[trans](statistic, options, extra_info)
        previous_stage=stage
    return statistic


def save_seqs(seqs, filenames='', seq_names=None):
    '''
    This gets a bit convoluted, because seqs is a list of different individuals.
    Each individual is itself a list of segments
    A segment is a list of basepairs.
    
    The output on the other hand is a different file for each segment
    Each file has a name-line and a sequence-line for each individual
    
    So the two formats are nested differently. The output format is saved in output_seqs
    '''
    
    try:
        basestring
    except NameError:
        basestring = str
    
    output_seqs=map(list, zip(*seqs))
    
    if isinstance(filenames, basestring):
        filenames=[filename+str(i+1) for i,filename in enumerate()]
    
    for filename,segment in zip(filenames,output_seqs):
        pass
        
        
    
def get_input(start_stage, options, extra_info={}):
    if start_stage==1:
        return get_populations(options.no_pops, options.population_names, extra_info)
    if start_stage==2:
        return read_ancestor_files(options.ancestor_files, extra_info)
    if start_stage==3:
        return read_allele_frequencies(options.allele_frequencies, options.truncate_af, extra_info)
    if start_stage==4:
        assert False, 'You should first read in the ancestor files before you can load from them'
    if start_stage==5:
        assert False, 'You should first read in the ancestor files before you can load from them'
    if start_stage==6:
        sequences_to_extract=get_wrapped_seqs(options.seq, options.seq_indices, options.ploidy_discrepancy)
        return read_sequences(options.seq_files, sequences_to_extract, extra_info)
    
def get_populations(no_pops, population_names, extra_info):
    if population_names:
        extra_info['pop_names']=population_names
        return population_names
    else:
        alphabet=list(map(chr, range(97, 123)))
        res=[]
        for i in range(no_pops):
            new_pop=alphabet[i]
            r=i//(123-97)
            if r>0:
                new_pop+=str(r)
            res.append(new_pop)
        extra_info['pop_names']=res
        return res
  
def get_wrapped_pops(pops, ancestor_multiplier, generations, pop_names):
    #assert len(pops)==ancestor_multiplier*2**generations, 'not enough true ancestors where specified'
    #check_list=list(set(pops))
    #assert all(pop in pop_names for pop in check_list), 'discrepancy between the true ancestors and the populations they could be chosen from'
    return [pops[x:x+2**generations] for x in range(0, len(pops), 2**generations)]

    

def get_wrapped_seqs(seq, seq_indices, ploidy_discrepancy):
    if len(seq)>0:
        s=seq
    elif len(seq_indices)>0:
        s=seq_indices
    else:
        assert seq, 'No sequences chosen'
    assert len(s)%ploidy_discrepancy==0, 'wrong number of sequences supplied'
    return [s[x:x+ploidy_discrepancy] for x in range(0, len(s), ploidy_discrepancy)]
    
def read_allele_frequencies(files, truncate_af, extra_info):
    all_allele_frequencies=[]
    setups=[]
    all_pop_names=[]
    for r in files:
        allele_frequencies={}
        lengths=[]
        pop_names=[]
        with open(r, "r") as f:
            text=list(f.readlines())
            for i,j in zip(text[0::2], text[1::2]):
                name=i.rstrip()
                freqs=list(map(float, j.split()))
                allele_frequencies[name]=[min(1-truncate_af, max(truncate_af,freq)) for freq in freqs ]
                pop_names.append(name)
            all_pop_names.append(pop_names)
            setups.append(len(freqs))
        all_allele_frequencies.append(allele_frequencies)
    if 'setups' in extra_info:
        print(setups)
        print(extra_info['setups'])
        assert setups==extra_info['setups'], 'the allele frequency file did not match some previously read files (recombination or ancestor or seq files)'
    else:
        extra_info['setups']=setups
    assert all(names==all_pop_names[0] for names in all_pop_names), 'some allele frequency files did not have identical populations.'
    extra_info['pop_names']=all_pop_names[0]
    extra_info['allele_frequencies']=all_allele_frequencies
    return all_allele_frequencies

def read_sequences(files, which_to_choose, extra_info):
    assert all(0<len(l)<=2 for l in which_to_choose), "One or two sequences are needed"
    all_seqs=[]
    all_names=[]
    setups=[]
    for r in files:
        haplotypes=[]
        names=[]
        with open(r, "r") as f:
            text=list(f.readlines())
            for i,j in zip(text[0::2], text[1::2]):
                name=i.rstrip()
                hap=list(map(int, j.split()))
                haplotypes.append(hap)
                names.append(name)
            assert all(len(h)==len(hap) for h in haplotypes), 'Some sequences were not equally long'
            setups.append(len(hap))
            all_names.append(names)
        all_seqs.append(haplotypes)
    if 'setups' in extra_info:
        assert setups==extra_info['setups'], 'the allele frequency file did not match some previously read files (recombination or ancestor)'
    else:
        extra_info['setups']=setups    
    assert all(names==all_names[0] for names in all_names), 'some ancestor files did not have identical ancestors.'
    extra_info['names']=all_names[0]
    names_to_indices={name:n for n,name in enumerate(all_names[0])}
#     print(len(all_names[0]), all_names[0])
#     for l in all_seqs:
#         print('\t',len(l))
#         for k in l:
#             print('\t'*2,len(k))
    
    chosen_seqs=[]
#     print('which_to_choose', which_to_choose)
    for segment_across_specimens in all_seqs:
        segment_across_samples=[]
        for t in which_to_choose:
            segment_within_sample=[]
            for r in t:
                if isinstance(r, int):
                    segment_within_sample.append(segment_across_specimens[r])
                else:
                    segment_within_sample.append(segment_across_specimens[names_to_indices[r]])
            if len(segment_within_sample)>1:
                segment_within_sample=[i+j for i,j in zip(segment_within_sample[0],segment_within_sample[1])]
            else:
                segment_within_sample=segment_within_sample[0]
            segment_across_samples.append(segment_within_sample)
        chosen_seqs.append(segment_across_samples)
#     print('chosen_seqs',chosen_seqs)
#     print('..............')
#     print(len(chosen_seqs))
#     for l in chosen_seqs:
#         print('\t',len(l))
#         for k in l:
#             print('\t'*2,len(k))
    #sorry about the readability in this one but it is just transposing:
    res=list(map(list,list(zip(*chosen_seqs))))
    print(len(res))
    for l in res:
        print('\t',len(l))
        for k in l:
            print('\t'*2,len(k))
    return res
    
    
def read_ancestor_files(files, extra_info):
    setups=[]
    all_ancestors=[]
    all_names=[]
    for r in files:
        haplotypes=[]
        names=[]
        with open(r, "r") as f:
            text=list(f.readlines())
            for i,j in zip(text[0::2], text[1::2]):
                name=i.rstrip()
                hap=list(map(int, j.split()))
                haplotypes.append(hap)
                names.append(name)
            assert all(len(h)==len(hap) for h in haplotypes), 'Some sequences were not equally long'
            setups.append(len(hap))
            all_names.append(names)
        all_ancestors.append(haplotypes)
    if 'setups' in extra_info:
        assert setups==extra_info['setups'], 'the recombination file and the ancestor file does not match in size'
    else:
        extra_info['setups']=setups
    assert all(names==all_names[0] for names in all_names), 'some ancestor files did not have identical ancestors.'
    #print('all_names',all_names[0])
    #print('files', files)
    extra_info['names']=all_names[0]
    return all_ancestors

if __name__=='__main__':
    from id_dic import id_dic
    class Object(object):
        pass

    options=Object()
    options.sequences_pipeline=[1,3,4,5,6,7,8]
    
    options.true_pops=None
    options.no_pops=2
    options.population_names=[]
    options.SNPs=37
    options.reps=1
    options.generations=3
    options.short_to_full=id_dic()
    options.simulate_recombination=True
    options.sim_recomb_rate=0.02
    options.recomb_rate=0.1
    options.skewness=1
    options.sequence_multiplier=1
    options.thin_coef=2
    ei={}
    
    a=get_seqs(options,ei)
    print(a)
    
    #print(get_wrapped_seqs(seq=[], seq_indices=[0,1], ploidy_discrepancy=1))
    
