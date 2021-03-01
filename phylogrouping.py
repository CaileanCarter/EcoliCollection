import numpy as np

def phylo(x):
    try:
        p = open(path.join(phylogroups, f"{x}.txt")).readlines()[-1]
        result = p.split()[1]
        return result
    except:
        return np.NaN

# phylo = lambda x : open(path.join(phylogroups, f"{x}.txt")).readlines()[-1]

ena["phylogroup"] = ena.index.map(phylo)