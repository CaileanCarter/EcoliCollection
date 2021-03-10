from os import path, system, mkdir
import glob
import argparse

import numpy as np
import pandas as pd

from fasta_from_ena import to_excel


def EzClerment(fp, output, script="EzClermont.py"):
    OUT = path.basename(fp).replace(".fasta", ".txt")
    system(f'python "{script}" "{fp}" 1> "{output}{OUT}" 2>&1')


def run_phylo(fp, script="EzClermont.py"):
    out = path.join(fp, "phylo")
    if not path.isdir(out):
        mkdir(out)
    for f in glob.glob(path.join(fp, "*.fasta")):
        EzClerment(fp=f, output=out, script=script)


def fetch_phylogroup(ena : pd.DataFrame, fp):

    def phylo(x):
        try:
            p = open(path.join(fp, f"{x}.txt")).readlines()[-1]
            result = p.split()[1]
            return result
        except FileNotFoundError:
            return np.NaN

    ena["phylogroup"] = ena.index.map(phylo)
    return ena


def fetch_serotype(ena : pd.DataFrame, fp):

    def ectype(x):
        PATH = path.join(fp, f"{x}.tabular")
        try:
            result = open(PATH).readlines()[-1].split()
        except FileNotFoundError:
            return np.NaN
        
        Otype = result[1]
        Htype = result[2]

        response = " ".join([Otype, Htype])

        return response

    ena["serotype"] = ena.index.map(ectype)
    return ena



def main():

    pass


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('-i', '--input', type=str)
    parser.add_argument('-s', '--serotype', type=str)
    parser.add_argument('-p', '--phylo', type=str)

    args = parser.parse_args()
    main()

    pass



if __name__ == "__main__":
    parse_arguments()


