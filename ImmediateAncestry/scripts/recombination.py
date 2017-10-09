
def read_recombination_file(files):
    lengths=[]
    recombs=[]
    setups=[]
    print(files)
    for r in files:
        with open(r, "r") as f:
            recombs.append(list(map(float, f.readline().split())))
        lengths.append(len(recombs[-1])+1)
    setups.append(lengths)
    return recombs, setups