def print_sequence_structure(sequences):
    print('no_of_individuals=',len(sequences))
    for i, seq_system in enumerate(sequences):
        print('individual number ',str(i+1), ' has ', len(seq_system), ' chromosomes of lengths')
        print('\t', map(str,map(len,seq_system)))
        
def print_recombination_structure(recombs):
    print('recombination structure, lenghts of chromosomes')
    print('\t', map(str,map(len,recombs)))
        
        