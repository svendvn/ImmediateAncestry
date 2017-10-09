def read_shortcuts(filename):
    full_to_short={}
    short_to_full={}
    with open(filename, 'r') as f:
        for l in f.readlines():
            fs=l.split()
            if len(fs)==2:
                full_to_short[fs[0]]=fs[1]
                short_to_full[fs[1]]=fs[0]
        
    return full_to_short, short_to_full