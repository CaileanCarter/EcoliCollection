import re
from itertools import product
from os import path, system

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

DIR = "."


def EzClerment(fp, output, script="EzClermont.py"):
    OUT = path.basename(fp).replace(".fasta", ".txt")
    system(f'python "{script}" "{fp}" 1> "{output}{OUT}" 2>&1')


def phylo(x):
    try:
        p = open(path.join(DIR, f"{x}.txt")).readlines()[-1]
        result = p.split()[1]
        return result
    except:
        return np.NaN


def fetch_phylogroup(ena : pd.DataFrame):
    ena["phylogroup"] = ena.index.map(phylo)
    return ena


def ectype(x):
    PATH = path.join(DIR, f"{x}.tabular")
    try:
        result = open(PATH).readlines()[-1].split()
    except:
        return np.NaN
    
    Otype = result[1]
    Htype = result[2]

    response = " ".join([Otype, Htype])

    return response


def fetch_serotype(ena : pd.DataFrame):
    ena["serotype"] = ena.index.map(ectype)
    return ena


#--------------------------------------------------
# plotting GPA from roary as a heatmap

class GPA:

    def __init__(self): raise TypeError("Object GPA cannot be initialised. Use GPA.plot()")

    def generate_clean_gene_set(self, GPA):
        GPA_filter = {}
        for row in GPA.itertuples():
            search = re.split(r"_\d+|\d+$", row[0])
            gene = search[0]
            if gene in GPA_filter.keys():
                GPA_filter[gene] = np.maximum(GPA_filter[gene], row[1:])
            else:
                GPA_filter[gene] = row[1:]
        return GPA_filter


    def similarity_score(self, colA, colB):
        result = np.where(colA == colB, True, False).sum()
        perc = (result / len(colA)) * 100
        return result, perc


    @classmethod
    def plot(cls, fp):
        gpa = pd.read_table(fp, index_col=0)
        tidy = cls.generate_clean_gene_set(cls, gpa)
        tidy_df = pd.DataFrame.from_dict(tidy, orient="index", columns=gpa.columns)

        df = pd.DataFrame(index=gpa.columns, columns=gpa.columns)
        df = df.fillna(value=0.0)
        for colA, colB in product(df.index, df.columns):
            if df.at[colB, colA] != 0.0:
                continue
            _, perc = cls.similarity_score(cls, tidy_df[colA], tidy_df[colB])
            df.at[colA, colB] = perc
            df.at[colB, colA] = perc

        sns.clustermap(df, cmap="mako")
        plt.show()