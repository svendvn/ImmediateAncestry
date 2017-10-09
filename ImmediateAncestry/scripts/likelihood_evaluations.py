

def evaluate_list(likelihood, configs, generations, short_to_full):
    configs=parse_configs(configs, generations, short_to_full)
    res_dic=[]
    for config in configs:
        print(config)
        res_dic.append((tuple(config),likelihood(config)))
    return res_dic
    
    
def parse_configs(unparsed_configs, generations, short_to_full):
    print('parsing', unparsed_configs)
    res=[]
    unfinished_config=[]
    for unparsed_config in unparsed_configs:
        if '.' in unparsed_config:
            with open(unparsed_config, 'r') as f:
                a=f.readlines()
                for l in a:
                    pieces=l.split()
                    if len(pieces)< 2**generations:
                        print('didnt read line', l)
                        continue
                    np=pieces[:2**generations]
                    print(np)
                    res.append([short_to_full[p] for p in np])
        else:
            pieces=unparsed_config.split()
            unfinished_config.extend(pieces)
            if len(unfinished_config)>=2**generations:
                np=unfinished_config[:2**generations]
                res.append([short_to_full[p] for p in np])
                unfinished_config=unfinished_config[2**generations:]
    return res

            