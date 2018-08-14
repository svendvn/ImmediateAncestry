from configuration_classes import find_smallest_equivalence_class
from argparse import ArgumentParser


def get_perm_number(candidate, true):
    if find_smallest_equivalence_class(candidate)==true:
        return 0
    for n1,s1 in enumerate(candidate):
        for n2, s2 in enumerate(candidate):
            candcand=list(candidate)
            candcand[n1]=s2
            candcand[n2]=s1
            if find_smallest_equivalence_class("".join(candcand))==true:
                if (n1<4 and n2<4) or (n1>=4 and n2>=4):
                    return 1
                else:
                    return 2
    return 3

def distance(candidate, true):
    max_num=len(candidate)
    if candidate==true:
        return "=",""
    hist1=dict((x, candidate.count(x)) for x in candidate)
    hist2=dict((x, true.count(x)) for x in true)
    no_diffs=0
    overscorer=""
    underscorer=""
    for v in set(list(hist1.keys())+list(hist2.keys())):
        k1=hist1.get(v,0)
        k2=hist2.get(v,0)
        if k1!=k2:
            if k1>k2:
                overscorer+=v
            else:
                underscorer+=v
            no_diffs+=abs(k1-k2)
    no_diffs=no_diffs//2
    ans=""
    if no_diffs==1:
        ans+="h1_"
        perm_numbers=[]
        for n in range(max_num):
            if candidate[n] in overscorer:
                candcand=list(candidate)
                candcand[n]=underscorer
                perm_numbers.append(get_perm_number("".join(candcand), true))
        perm_number=min(perm_numbers)
    elif no_diffs==0:
        perm_number=get_perm_number(candidate, true)
    elif no_diffs>=2:
        return 'h'+str(no_diffs)+'?', overscorer
    converter=["=","1w","1b","+"]
    ans+=converter[perm_number]
    return ans,overscorer

def get_config_from_file(config_file):
    with open(config_file,'r') as f:
        a=f.readline().rstrip().split()
        subspecies=''.join(a[:-1])
    return subspecies
    
def get_config_from_filename_aug13(config_file):
    return config_file.split('_')[2]

def get_info_from_file_name_aug13(config_file):
    
    res={}
    parts=config_file.split('_')
    res['complexity']=int(parts[0][1:])
    res['true_config']=parts[1]
    res['generations']=int(parts[3][1:])
    res['bin_size']=int(parts[4][1:])
    res['simulation']=int(parts[5][1:])
    res['id_number']=int(parts[6][1:].split('.')[0])
    return res

def write_to_file(dic, filename='tmp.out'):
    ad=sorted(list(dic.items()))
    with open(filename,'w') as f:
        f.write(','.join((str(a[1]) for a in ad)))
    

if __name__=='__main__':
    
    usage='parse simulation study from aug13, 2018.'

    
    parser = ArgumentParser(usage=usage)
    
    #input/output arguments
    parser.add_argument("--config_file", type=str, nargs='*', default=['s0_tttt_0_g2_b10_w3_i55.out'], help='file_to_read_From')
    parser.add_argument('--output_file_suffix', type=str, default='2', help='file to write results to') 
                        
    options=parser.parse_args()
    
    for conf_file in options.config_file:
        info_dic=get_info_from_file_name_aug13(conf_file)
        info_dic['inferred']=get_config_from_file(conf_file)
        dist=distance(find_smallest_equivalence_class(info_dic['inferred']),
                      find_smallest_equivalence_class(info_dic['true_config']))
        info_dic['distance']=str(dist[0])
        info_dic['overrepresented']=str(dist[1])
        write_to_file(info_dic, conf_file+options.output_file_suffix)
        
    print(sorted(info_dic.keys()))

    
    print(distance("eeteeeee","eeeeeeee"))