

def evaluate_list(likelihood, configs, generations, short_to_full):
    print ("#evaluate_list(likelihood, configs, generations, short_to_full)")
    configs=parse_configs(configs, generations, short_to_full)
    res_dic=[]
    for config in configs:
        print("config = ", config)
        res_dic.append((tuple(config),likelihood(config)))
    print ("return res_dic", res_dic,"\n")
    return res_dic
    
    
def parse_configs(unparsed_configs, generations, short_to_full):
    print ( "#parse_configs")
    print('unparsed_configs', unparsed_configs)
    res=[]
    unfinished_config=[]
    for unparsed_config in unparsed_configs:
        if '.' in unparsed_config:
            print ("'.' in unparsed_config")
            with open(unparsed_config, 'r') as f:
                a=f.readlines()
                print ("a = " ,a)
                for l in a:
                    pieces=l.split()
                    print ("pieces = ", pieces)
                    if len(pieces)< 2**generations:
                        print('didnt read line', l)
                        continue
                    np=pieces[:2**generations]
                    print("np = ", np)
                    res.append([short_to_full[p] for p in np])
        else:
            print ("'.' No in unparsed_config")
            pieces=unparsed_config.split()
            print ("pieces = ", pieces)
            unfinished_config.extend(pieces)
            if len(unfinished_config)>=2**generations:
                print ( "len(unfinished_config)>=2**generations")
                np=unfinished_config[:2**generations]
                print ("np = ", np)
                res.append([short_to_full[p] for p in np])
                unfinished_config=unfinished_config[2**generations:]
                print ("unfinished_config", unfinished_config)
    return res

            
