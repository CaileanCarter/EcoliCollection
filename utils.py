from os import path, system

import numpy as np
import pandas as pd

DIR = "."


def EzClerment(fp, output, script="EzClermont.py"):
    OUT = path.basename(fp).replace(".fasta", ".txt")
    system(f'python "{script}" "{fp}" 1> "{output}{OUT}" 2>&1')


def phylo(x):
    try:
        p = open(path.join(DIR, f"{x}.txt")).readlines()[-1]
        result = p.split()[1]
        return result
    except FileNotFoundError:
        return np.NaN


def fetch_phylogroup(ena : pd.DataFrame):
    ena["phylogroup"] = ena.index.map(phylo)
    return ena


def ectype(x):
    PATH = path.join(DIR, f"{x}.tabular")
    try:
        result = open(PATH).readlines()[-1].split()
    except FileNotFoundError:
        return np.NaN
    
    Otype = result[1]
    Htype = result[2]

    response = " ".join([Otype, Htype])

    return response


def fetch_serotype(ena : pd.DataFrame):
    ena["serotype"] = ena.index.map(ectype)
    return ena


